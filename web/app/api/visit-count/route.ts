import { NextResponse } from "next/server";

export const runtime = "nodejs";

type UpstashNumberResponse = {
  result?: number | string | null;
};

declare global {
  var nifFinderLocalVisitCount: number | undefined;
}

const visitCountKey = "nif-finder:visit-count";

async function upstashCommand(command: "get" | "incr") {
  const redisUrl = process.env.UPSTASH_REDIS_REST_URL;
  const redisToken = process.env.UPSTASH_REDIS_REST_TOKEN;
  if (!redisUrl || !redisToken) return null;

  const response = await fetch(`${redisUrl}/${command}/${visitCountKey}`, {
    headers: { authorization: `Bearer ${redisToken}` },
    cache: "no-store",
  });
  if (!response.ok) {
    throw new Error(`Visit counter storage returned ${response.status}.`);
  }

  const payload = (await response.json()) as UpstashNumberResponse;
  const count = Number(payload.result ?? 0);
  return Number.isFinite(count) ? count : 0;
}

function localCount(increment: boolean) {
  if (process.env.NODE_ENV === "production") return null;
  globalThis.nifFinderLocalVisitCount ??= 0;
  if (increment) {
    globalThis.nifFinderLocalVisitCount += 1;
  }
  return globalThis.nifFinderLocalVisitCount;
}

async function visitCount(increment: boolean) {
  try {
    const upstashCount = await upstashCommand(increment ? "incr" : "get");
    if (upstashCount !== null) {
      return NextResponse.json({ enabled: true, count: upstashCount, provider: "upstash" });
    }

    const count = localCount(increment);
    if (count !== null) {
      return NextResponse.json({ enabled: true, count, provider: "local" });
    }

    return NextResponse.json({ enabled: false });
  } catch (error) {
    return NextResponse.json(
      {
        enabled: false,
        error: error instanceof Error ? error.message : String(error),
      },
      { status: 200 },
    );
  }
}

export async function GET() {
  return visitCount(false);
}

export async function POST() {
  return visitCount(true);
}
