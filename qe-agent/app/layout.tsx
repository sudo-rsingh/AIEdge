import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "QE Input Generator - AI Agent",
  description: "AI-powered Quantum Espresso input file generator using LangChain",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body className="antialiased">
        {children}
      </body>
    </html>
  );
}
