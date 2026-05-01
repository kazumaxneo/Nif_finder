import { NextRequest, NextResponse } from "next/server";
import { randomUUID } from "node:crypto";
import { mkdir, readFile, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import path from "node:path";
import { promisify } from "node:util";
import { execFile } from "node:child_process";

export const runtime = "nodejs";
export const maxDuration = 60;

const execFileAsync = promisify(execFile);

type AnalyzeRequest = {
  fasta?: string;
  jobs?: number;
  cpu?: number;
  plot?: boolean;
  evalue?: number;
};

function countFastaRecords(fasta: string) {
  return fasta
    .split(/\r?\n/)
    .filter((line) => line.trim().startsWith(">")).length;
}

function parseNifFinderTsv(tsv: string) {
  return tsv
    .trim()
    .split(/\r?\n/)
    .slice(1)
    .filter(Boolean)
    .map((line) => {
      const [query, logEvalue, alignLength, queryLength, prediction, completeness] = line.split("\t");
      return {
        query,
        logEvalue: Number(logEvalue),
        alignLength: Number(alignLength),
        queryLength: Number(queryLength),
        prediction,
        completeness,
      };
    });
}

async function runLocalNifFinder(fasta: string, jobs: number, cpu: number, plot: boolean, evalue: number) {
  const root = process.env.NIF_FINDER_ROOT
    ? path.resolve(process.env.NIF_FINDER_ROOT)
    : path.resolve(process.cwd(), "..", "generl_bacteria");
  const python = process.env.NIF_FINDER_PYTHON ?? "python";
  const script = path.join(root, "Nif_finderv0_24.py");
  const workDir = path.join(tmpdir(), `nif-finder-web-${randomUUID()}`);
  const queryPath = path.join(workDir, "query.faa");
  const outPrefix = path.join(workDir, "result");

  await mkdir(workDir, { recursive: true });
  try {
    await writeFile(queryPath, fasta);
    const args = [
      script,
      "-q",
      queryPath,
      "-o",
      outPrefix,
      "--jobs",
      String(jobs),
      "--cpu",
      String(cpu),
      "--evalue",
      String(evalue),
    ];
    if (plot) {
      args.push("-p");
    }

    await execFileAsync(
      python,
      args,
      {
        env: {
          ...process.env,
          NIF_FINDER_DB: process.env.NIF_FINDER_DB ?? root,
        },
        timeout: 55_000,
      },
    );
    const tsv = await readFile(`${outPrefix}.txt`, "utf8");
    let plotPngBase64: string | null = null;
    try {
      const plot = await readFile(`${outPrefix}_scatter.png`);
      plotPngBase64 = plot.toString("base64");
    } catch {
      plotPngBase64 = null;
    }
    return { records: parseNifFinderTsv(tsv), plotPngBase64, runner: "local" };
  } finally {
    await rm(workDir, { recursive: true, force: true });
  }
}

export async function POST(request: NextRequest) {
  let body: AnalyzeRequest;
  try {
    body = (await request.json()) as AnalyzeRequest;
  } catch {
    return NextResponse.json({ error: "Invalid JSON request." }, { status: 400 });
  }

  const fasta = body.fasta?.trim();
  if (!fasta) {
    return NextResponse.json({ error: "Please provide a protein FASTA sequence." }, { status: 400 });
  }
  if (!fasta.startsWith(">")) {
    return NextResponse.json({ error: "Input must be FASTA formatted and start with '>'." }, { status: 400 });
  }

  const apiUrl = process.env.NIF_FINDER_API_URL;
  if (!apiUrl) {
    if (process.env.NIF_FINDER_LOCAL_RUNNER === "1") {
      try {
        const result = await runLocalNifFinder(
          fasta,
          body.jobs ?? 3,
          body.cpu ?? 6,
          body.plot ?? true,
          body.evalue ?? 1e-10,
        );
        return NextResponse.json(result);
      } catch (error) {
        return NextResponse.json(
          {
            error: "Local Nif-Finder runner failed.",
            detail: error instanceof Error ? error.message : String(error),
          },
          { status: 500 },
        );
      }
    }

    return NextResponse.json(
      {
        error: "NIF_FINDER_API_URL is not configured for this deployment.",
        setup:
          "Deploy or expose a compute service that runs the Python/HMMER Nif-Finder pipeline, then set NIF_FINDER_API_URL in Vercel. For local development, set NIF_FINDER_LOCAL_RUNNER=1.",
        sequenceCount: countFastaRecords(fasta),
      },
      { status: 501 },
    );
  }

  let upstream: Response;
  const headers: Record<string, string> = { "content-type": "application/json" };
  if (process.env.NIF_FINDER_API_KEY) {
    headers["x-api-key"] = process.env.NIF_FINDER_API_KEY;
  }

  try {
    upstream = await fetch(apiUrl, {
      method: "POST",
      headers,
      body: JSON.stringify({
        fasta,
        jobs: body.jobs ?? 3,
        cpu: body.cpu ?? 6,
        plot: body.plot ?? true,
        evalue: body.evalue ?? 1e-10,
      }),
    });
  } catch (error) {
    return NextResponse.json(
      {
        error: "Could not reach the Nif-Finder compute API.",
        detail: error instanceof Error ? error.message : String(error),
      },
      { status: 502 },
    );
  }

  const text = await upstream.text();
  let payload: unknown;
  try {
    payload = text ? JSON.parse(text) : { error: "Compute API returned an empty response." };
  } catch {
    payload = {
      error: "Compute API returned a non-JSON response.",
      detail: text.slice(0, 1000),
    };
  }

  return NextResponse.json(payload, { status: upstream.status });
}
