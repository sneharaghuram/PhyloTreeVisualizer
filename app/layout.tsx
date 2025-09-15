import "./globals.css";
import Header from "@/components/Header";
import Footer from "@/components/Footer";

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body className="bg-gray-100 flex flex-col min-h-screen">
        <Header />
        <main className="text-black flex-grow pt-13">{children}</main>
        <Footer />
      </body>
    </html>
  );
}
