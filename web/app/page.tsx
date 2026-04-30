"use client";

import { ChangeEvent, useMemo, useState } from "react";
import { Activity, AlertCircle, BarChart3, FileUp, Play, Settings } from "lucide-react";
import {
  CartesianGrid,
  Legend,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

type ResultRecord = {
  query: string;
  logEvalue: number;
  alignLength: number;
  queryLength: number;
  prediction: string;
  completeness: string;
};

type ApiResponse = {
  records?: ResultRecord[];
  error?: string;
  setup?: string;
  sequenceCount?: number;
};

const sampleFasta = `>Calothrix_sp.NIES-4101_01433
MAVKIAINGFGRIGRNVRAAQKRLKAEGKVVLVTGKGGIGKSTTSQNTLAALGVKVLQIGCDPKHDSTFTLTGHGEGAPEGTVKALGEAVGLGAEVVKV
>Calothrix_sp.NIES-4101_03559
MSIVNIGPGQIGKSTTSQNTLAALAQGLNVLNVGCDPKHDSTFTLTGAGQGAPEGTVKAIGEAVGLGVEVVKVTE
`;

const geneColors: Record<string, string> = {
  nifH: "#d73a31",
  nifD: "#206fb1",
  nifK: "#1f8f55",
  nifE: "#b56b00",
  nifN: "#7b4fb3",
  nifB: "#0f8b8d",
  other: "#7f8790",
  unclassifiable: "#9aa0a6",
};

export default function Home() {
  const [fasta, setFasta] = useState(sampleFasta);
  const [jobs, setJobs] = useState(3);
  const [cpu, setCpu] = useState(6);
  const [loading, setLoading] = useState(false);
  const [response, setResponse] = useState<ApiResponse | null>(null);

  const records = response?.records ?? [];
  const geneCounts = useMemo(() => {
    return records.reduce<Record<string, number>>((acc, record) => {
      acc[record.prediction] = (acc[record.prediction] ?? 0) + 1;
      return acc;
    }, {});
  }, [records]);

  async function analyze() {
    setLoading(true);
    setResponse(null);
    try {
      const res = await fetch("/api/analyze", {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ fasta, jobs, cpu }),
      });
      const data = (await res.json()) as ApiResponse;
      setResponse(data);
    } catch (error) {
      setResponse({ error: error instanceof Error ? error.message : "Request failed." });
    } finally {
      setLoading(false);
    }
  }

  function handleFile(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => setFasta(String(reader.result ?? ""));
    reader.readAsText(file);
  }

  return (
    <main className="workspace">
      <aside className="sidebar">
        <div className="brand">
          <Activity size={26} aria-hidden />
          <div>
            <h1>Nif-Finder Web</h1>
            <p>Protein FASTA analysis and nif hit visualization</p>
          </div>
        </div>

        <label className="field">
          Protein FASTA
          <textarea value={fasta} onChange={(event) => setFasta(event.target.value)} spellCheck={false} />
        </label>

        <div className="upload-row">
          <label className="file-button" title="Load FASTA from a local file">
            <FileUp size={18} aria-hidden />
            <span>Load FASTA</span>
            <input type="file" accept=".faa,.fa,.fasta,.txt" onChange={handleFile} />
          </label>
          <button className="ghost-button" type="button" onClick={() => setFasta(sampleFasta)}>
            Sample
          </button>
        </div>

        <div className="settings">
          <div className="section-title">
            <Settings size={17} aria-hidden />
            Run settings
          </div>
          <label>
            Jobs
            <input type="number" min={1} max={6} value={jobs} onChange={(event) => setJobs(Number(event.target.value))} />
          </label>
          <label>
            CPU
            <input type="number" min={1} max={24} value={cpu} onChange={(event) => setCpu(Number(event.target.value))} />
          </label>
        </div>

        <button className="run-button" type="button" onClick={analyze} disabled={loading}>
          <Play size={18} aria-hidden />
          {loading ? "Running" : "Run analysis"}
        </button>
      </aside>

      <section className="results">
        <div className="summary-row">
          <div>
            <p className="eyebrow">Results</p>
            <h2>Nif hit overview</h2>
          </div>
          <div className="counter">{records.length} hits</div>
        </div>

        {response?.error ? (
          <div className="notice">
            <AlertCircle size={20} aria-hidden />
            <div>
              <strong>{response.error}</strong>
              {response.setup ? <p>{response.setup}</p> : null}
              {response.sequenceCount ? <p>{response.sequenceCount} FASTA record(s) were parsed.</p> : null}
            </div>
          </div>
        ) : null}

        {records.length > 0 ? (
          <>
            <div className="metrics">
              {Object.entries(geneCounts).map(([gene, count]) => (
                <div className="metric" key={gene}>
                  <span style={{ backgroundColor: geneColors[gene] ?? geneColors.other }} />
                  <div>
                    <strong>{count}</strong>
                    <p>{gene}</p>
                  </div>
                </div>
              ))}
            </div>

            <div className="chart-panel">
              <div className="section-title">
                <BarChart3 size={18} aria-hidden />
                Protein length vs -log10(E-value)
              </div>
              <ResponsiveContainer width="100%" height={360}>
                <ScatterChart margin={{ top: 16, right: 16, bottom: 24, left: 4 }}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis dataKey="queryLength" name="Protein length" unit=" aa" />
                  <YAxis dataKey="logEvalue" name="-log10(E)" />
                  <Tooltip cursor={{ strokeDasharray: "3 3" }} />
                  <Legend />
                  {Object.keys(geneCounts).map((gene) => (
                    <Scatter
                      key={gene}
                      name={gene}
                      data={records.filter((record) => record.prediction === gene)}
                      fill={geneColors[gene] ?? geneColors.other}
                    />
                  ))}
                </ScatterChart>
              </ResponsiveContainer>
            </div>

            <div className="table-wrap">
              <table>
                <thead>
                  <tr>
                    <th>Query</th>
                    <th>Prediction</th>
                    <th>Completeness</th>
                    <th>-log E</th>
                    <th>Align len</th>
                    <th>Protein len</th>
                  </tr>
                </thead>
                <tbody>
                  {records.map((record) => (
                    <tr key={`${record.query}-${record.prediction}`}>
                      <td>{record.query}</td>
                      <td>{record.prediction}</td>
                      <td>{record.completeness}</td>
                      <td>{record.logEvalue.toFixed(2)}</td>
                      <td>{record.alignLength}</td>
                      <td>{record.queryLength}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </>
        ) : (
          <div className="empty-state">
            Submit a protein FASTA file to view predicted nif hits, completeness calls, and scatter plot output.
          </div>
        )}
      </section>
    </main>
  );
}
