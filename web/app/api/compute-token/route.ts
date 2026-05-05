import { createHmac, randomUUID } from "node:crypto";
import { NextResponse } from "next/server";

export const runtime = "nodejs";

export async function POST() {
  const apiUrl = process.env.NIF_FINDER_API_URL;
  const apiKey = process.env.NIF_FINDER_API_KEY;
  if (!apiUrl || !apiKey) {
    return NextResponse.json({ error: "Compute API token signing is not configured." }, { status: 501 });
  }

  const expires = String(Math.floor(Date.now() / 1000) + 600);
  const nonce = randomUUID();
  const signature = createHmac("sha256", apiKey).update(`${expires}.${nonce}`).digest("hex");

  return NextResponse.json({
    apiUrl,
    token: `${expires}.${nonce}.${signature}`,
  });
}
