import "dotenv/config";
import express from "express";
import path from "path";

async function startServer() {
  const app = express();
  const PORT = 3000;

  app.use(express.json({ limit: '50mb' }));

  // API routes FIRST
  app.get("/api/health", (req, res) => {
    res.json({ status: "ok" });
  });

  app.post("/api/diffdock", async (req, res) => {
    try {
      const apiKey = process.env.NVIDIA_API_KEY;
      if (!apiKey) {
        return res.status(500).json({ error: "NVIDIA_API_KEY is not configured" });
      }

      console.log("Downloading protein from RCSB...");
      const proteinResponse = await fetch("https://files.rcsb.org/download/8G43.pdb");
      if (!proteinResponse.ok) throw new Error("Failed to download protein");
      const proteinText = await proteinResponse.text();
      const proteinContent = proteinText.split("\n").filter(line => line.startsWith("ATOM")).join("\n");

      console.log("Downloading ligand from RCSB...");
      const ligandResponse = await fetch("https://files.rcsb.org/ligands/download/ZU6_ideal.sdf");
      if (!ligandResponse.ok) throw new Error("Failed to download ligand");
      const ligandContent = await ligandResponse.text();

      const url = "https://health.api.nvidia.com/v1/biology/mit/diffdock";
      
      console.log("Making request to DiffDock...");
      const response = await fetch(url, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "Authorization": `Bearer ${apiKey}`
        },
        body: JSON.stringify({
          ligand: ligandContent,
          ligand_file_type: "sdf",
          protein: proteinContent,
          num_poses: 10,
          time_divisions: 20,
          steps: 18,
          save_trajectory: false
        })
      });

      if (!response.ok) {
        const errorText = await response.text();
        console.error("DiffDock API error:", errorText);
        return res.status(response.status).json({ error: `DiffDock API error: ${response.statusText}` });
      }

      const data = await response.json();
      res.json(data);
    } catch (error) {
      console.error("Error calling DiffDock:", error);
      res.status(500).json({ error: "Internal server error" });
    }
  });

  // Vite middleware for development
  if (process.env.NODE_ENV !== "production") {
    const { createServer: createViteServer } = await import("vite");
    const vite = await createViteServer({
      server: { middlewareMode: true },
      appType: "spa",
    });
    app.use(vite.middlewares);
  } else {
    const distPath = path.join(process.cwd(), 'dist');
    app.use(express.static(distPath));
    app.get('*all', (req, res) => {
      res.sendFile(path.join(distPath, 'index.html'));
    });
  }

  app.listen(PORT, "0.0.0.0", () => {
    console.log(`Server running on http://localhost:${PORT}`);
  });
}

startServer();
