import { NextRequest, NextResponse } from "next/server";

export const runtime = "nodejs";
export const maxDuration = 60;

function computeEndpoint(pathname: string) {
  const apiUrl = process.env.NIF_FINDER_API_URL;
  if (!apiUrl) return null;
  const url = new URL(apiUrl);
  url.pathname = url.pathname.replace(/\/analyze\/?$/, `/${pathname}`).replace(/\/$/, `/${pathname}`);
  if (!url.pathname.endsWith(`/${pathname}`)) {
    url.pathname = `/${pathname}`;
  }
  return url.toString();
}

export async function POST(request: NextRequest) {
  const endpoint = computeEndpoint("compare-clusters");
  if (!endpoint) {
    return NextResponse.json(
      { error: "NIF_FINDER_API_URL is not configured for this deployment." },
      { status: 501 },
    );
  }

  let body: unknown;
  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON request." }, { status: 400 });
  }

  const headers: Record<string, string> = { "content-type": "application/json" };
  if (process.env.NIF_FINDER_API_KEY) {
    headers["x-api-key"] = process.env.NIF_FINDER_API_KEY;
  }

  let upstream: Response;
  try {
    upstream = await fetch(endpoint, {
      method: "POST",
      headers,
      body: JSON.stringify(body),
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
