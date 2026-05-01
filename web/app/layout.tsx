import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Nif-Finder Web",
  description: "Submit protein FASTA sequences and visualize Nif-Finder results.",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
