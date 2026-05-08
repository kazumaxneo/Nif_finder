"use client";

import { ChangeEvent, useMemo, useState } from "react";
import { AlertCircle, Download, FileUp, Play, Settings } from "lucide-react";

type ResultRecord = {
  query: string;
  logEvalue: number;
  alignLength: number;
  queryLength: number;
  prediction: string;
  completeness: string;
  operonLabel?: string | null;
};

type ApiResponse = {
  records?: ResultRecord[];
  plotPngBase64?: string | null;
  genomicContextOverviewSvg?: string | null;
  genomicContextLocalSvg?: string | null;
  genomicContextMessage?: string | null;
  error?: string;
  setup?: string;
  detail?: string;
  sequenceCount?: number;
};

const nifGenes = ["nifH", "nifD", "nifK", "nifE", "nifN", "nifB"];
const exampleDatasets = [
  { id: "none", label: "None" },
  { id: "leptolyngbya-boryana-dg5", label: "Leptolyngbya boryana dg5" },
  { id: "calothrix-fragmented", label: "Calothrix strain's fragmented nif genes" },
];

const figure1Caption = (
  <>
    Figure 1. 2D Similarity Plot of homology search for the six Nif proteins encoded by <em>nifHDKENB</em> of
    various bacterial genomes. Background reference points show the relationship between –log10(E-value) from HMMER3
    searches using <em>nif</em> HMM profiles and hit protein length across proteomes of 586 cyanobacterial strains and
    other bacterial strains. Query sequences submitted to Nif-Finder are overlaid and highlighted on this 2D Similarity
    Plot, and are classified by their best matching <em>nif</em> profile. The color indicates the single best hit
    protein from SWISS-PROT sequences. Circle, triangle, and star plots represent hits from complete genome assemblies,
    draft genome assemblies, and Nif fusion proteins, respectively.
  </>
);

const figure2Caption = (
  <>
    Figure 2. Genomic locations of <em>nif</em> genes identified by Nif-Finder. Colored arrows indicate <em>nif</em>{" "}
    genes identified by Nif-Finder, and gray arrows indicate neighboring coding sequences. Enlarged regional views show
    local genomic neighborhoods around <em>nif</em> hits. A whole-genome view is shown when it provides additional
    positional context and is omitted when redundant. Labels below highlighted arrows indicate hit status: full-length,
    fragment, or operon.
  </>
);

const table2Caption = (
  <>
    Table 2. Detailed list of Nif-Finder hits. Query indicates the input protein identifier, Prediction indicates the
    best Nif-Finder assignment, Completeness indicates full-length, operon, fragment, or unclassifiable status, and the
    remaining columns report HMMER search score and protein length metrics.
  </>
);

type FastaEntry = {
  id: string;
  header: string;
  sequence: string;
};

type ZipEntry = {
  name: string;
  data: Uint8Array;
};

export default function Home() {
  const [fasta, setFasta] = useState("");
  const [genbank, setGenbank] = useState("");
  const [genbankFileName, setGenbankFileName] = useState("");
  const [jobs, setJobs] = useState(3);
  const [cpu, setCpu] = useState(6);
  const [plotOutput, setPlotOutput] = useState(true);
  const [showOnlyNifHits, setShowOnlyNifHits] = useState(false);
  const [exampleDataset, setExampleDataset] = useState("none");
  const [evalueThreshold, setEvalueThreshold] = useState("1e-10");
  const [loading, setLoading] = useState(false);
  const [response, setResponse] = useState<ApiResponse | null>(null);

  const records = response?.records ?? [];
  const displayedRecords = showOnlyNifHits
    ? records.filter((record) => nifGenes.includes(record.prediction))
    : records;
  const nifSummary = useMemo(() => {
    return nifGenes.map((gene) => {
      const geneRecords = records.filter((record) => record.prediction === gene);
      const full = geneRecords.filter((record) => record.completeness === "Full").length;
      const operon = geneRecords.filter((record) => record.completeness === "Full_operon").length;
      const incomplete = geneRecords.filter((record) => record.completeness === "Fragment").length;

      return {
        gene,
        total: geneRecords.length,
        full,
        operon,
        incomplete,
      };
    });
  }, [records]);
  const totalNifCopies = nifSummary.reduce((sum, row) => sum + row.total, 0);

  async function loadExampleDataset(dataset: string) {
    setExampleDataset(dataset);
    setResponse(null);
    if (dataset === "leptolyngbya-boryana-dg5") {
      const [fastaResponse, genbankResponse] = await Promise.all([
        fetch("/examples/leptolyngbya_boryana_dg5.faa"),
        fetch("/examples/leptolyngbya_boryana_dg5.gbk"),
      ]);
      setFasta(await fastaResponse.text());
      setGenbank(await genbankResponse.text());
      setGenbankFileName("leptolyngbya_boryana_dg5.gbk");
    } else if (dataset === "calothrix-fragmented") {
      const exampleResponse = await fetch("/examples/calothrix_fragmented_nif_genes.faa");
      setFasta(await exampleResponse.text());
      setGenbank("");
      setGenbankFileName("");
    } else {
      setFasta("");
      setGenbank("");
      setGenbankFileName("");
    }
  }

  function downloadText(filename: string, content: string, type = "text/plain;charset=utf-8") {
    const blob = new Blob([content], { type });
    downloadBlob(filename, blob);
  }

  function downloadBlob(filename: string, blob: Blob) {
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    link.click();
    URL.revokeObjectURL(url);
  }

  function resultRows(selectedRecords = displayedRecords) {
    return selectedRecords.map((record) => [
      record.query,
      record.logEvalue.toFixed(2),
      String(record.alignLength),
      record.queryLength == null ? "N/A" : String(record.queryLength),
      record.prediction,
      record.completeness,
      record.operonLabel ?? "",
    ]);
  }

  function resultsTsv(selectedRecords = displayedRecords) {
    const header = ["Query", "-log_Evalue", "Align_Len", "Query_Length", "Prediction", "Completeness", "Operon"];
    return [header, ...resultRows(selectedRecords)].map((row) => row.join("\t")).join("\n") + "\n";
  }

  function csvCell(value: string) {
    return /[",\n\r]/.test(value) ? `"${value.replace(/"/g, '""')}"` : value;
  }

  function resultsCsv(selectedRecords = displayedRecords) {
    const header = ["Query", "-log_Evalue", "Align_Len", "Query_Length", "Prediction", "Completeness", "Operon"];
    return [header, ...resultRows(selectedRecords)].map((row) => row.map(csvCell).join(",")).join("\n") + "\n";
  }

  function parseFastaEntries(input: string) {
    const entries: FastaEntry[] = [];
    let header = "";
    let sequenceParts: string[] = [];

    for (const rawLine of input.split(/\r?\n/)) {
      const line = rawLine.trim();
      if (!line) continue;
      if (line.startsWith(">")) {
        if (header) {
          entries.push({
            id: header.split(/\s+/)[0],
            header,
            sequence: sequenceParts.join(""),
          });
        }
        header = line.slice(1).trim();
        sequenceParts = [];
      } else if (header) {
        sequenceParts.push(line.replace(/\s+/g, ""));
      }
    }

    if (header) {
      entries.push({
        id: header.split(/\s+/)[0],
        header,
        sequence: sequenceParts.join(""),
      });
    }

    return entries;
  }

  function formatFasta(entries: FastaEntry[], selectedRecords: Map<string, ResultRecord>) {
    return entries
      .map((entry) => {
        const record = selectedRecords.get(entry.id);
        const annotation = record ? ` |Nif-Finder ${displayPrediction(record)} ${displayCompleteness(record)}` : "";
        const wrappedSequence = entry.sequence.match(/.{1,80}/g)?.join("\n") ?? "";
        return `>${entry.header}${annotation}\n${wrappedSequence}`;
      })
      .join("\n");
  }

  function fastaContent(selectedRecords: ResultRecord[]) {
    const recordMap = new Map<string, ResultRecord>();
    for (const record of selectedRecords) {
      if (!recordMap.has(record.query)) {
        recordMap.set(record.query, record);
      }
    }

    const entries = parseFastaEntries(fasta).filter((entry) => recordMap.has(entry.id));
    return entries.length === 0 ? "" : `${formatFasta(entries, recordMap)}\n`;
  }

  function displayPrediction(record: ResultRecord) {
    return record.completeness === "Full_operon" && record.operonLabel ? record.operonLabel : record.prediction;
  }

  function displayCompleteness(record: ResultRecord) {
    return record.completeness === "Full_operon" ? "Operon" : record.completeness;
  }

  function textBytes(text: string) {
    return new TextEncoder().encode(text);
  }

  function base64Bytes(base64: string) {
    const binary = atob(base64);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i += 1) {
      bytes[i] = binary.charCodeAt(i);
    }
    return bytes;
  }

  function crc32(data: Uint8Array) {
    let crc = 0xffffffff;
    for (const byte of data) {
      crc ^= byte;
      for (let i = 0; i < 8; i += 1) {
        crc = crc & 1 ? (crc >>> 1) ^ 0xedb88320 : crc >>> 1;
      }
    }
    return (crc ^ 0xffffffff) >>> 0;
  }

  function dosDateTime(date = new Date()) {
    const time = (date.getHours() << 11) | (date.getMinutes() << 5) | Math.floor(date.getSeconds() / 2);
    const day = Math.max(1, date.getDate());
    const dosDate = ((date.getFullYear() - 1980) << 9) | ((date.getMonth() + 1) << 5) | day;
    return { time, date: dosDate };
  }

  function concatBytes(parts: Uint8Array[]) {
    const totalLength = parts.reduce((sum, part) => sum + part.length, 0);
    const output = new Uint8Array(totalLength);
    let offset = 0;
    for (const part of parts) {
      output.set(part, offset);
      offset += part.length;
    }
    return output;
  }

  function zipHeader(size: number, write: (view: DataView) => void) {
    const bytes = new Uint8Array(size);
    const view = new DataView(bytes.buffer);
    write(view);
    return bytes;
  }

  function createZip(entries: ZipEntry[]) {
    const fileParts: Uint8Array[] = [];
    const centralParts: Uint8Array[] = [];
    let offset = 0;
    const { time, date } = dosDateTime();

    for (const entry of entries) {
      const nameBytes = textBytes(entry.name);
      const checksum = crc32(entry.data);
      const localHeader = zipHeader(30, (view) => {
        view.setUint32(0, 0x04034b50, true);
        view.setUint16(4, 20, true);
        view.setUint16(6, 0, true);
        view.setUint16(8, 0, true);
        view.setUint16(10, time, true);
        view.setUint16(12, date, true);
        view.setUint32(14, checksum, true);
        view.setUint32(18, entry.data.length, true);
        view.setUint32(22, entry.data.length, true);
        view.setUint16(26, nameBytes.length, true);
        view.setUint16(28, 0, true);
      });
      fileParts.push(localHeader, nameBytes, entry.data);

      const centralHeader = zipHeader(46, (view) => {
        view.setUint32(0, 0x02014b50, true);
        view.setUint16(4, 20, true);
        view.setUint16(6, 20, true);
        view.setUint16(8, 0, true);
        view.setUint16(10, 0, true);
        view.setUint16(12, time, true);
        view.setUint16(14, date, true);
        view.setUint32(16, checksum, true);
        view.setUint32(20, entry.data.length, true);
        view.setUint32(24, entry.data.length, true);
        view.setUint16(28, nameBytes.length, true);
        view.setUint16(30, 0, true);
        view.setUint16(32, 0, true);
        view.setUint16(34, 0, true);
        view.setUint16(36, 0, true);
        view.setUint32(38, 0, true);
        view.setUint32(42, offset, true);
      });
      centralParts.push(centralHeader, nameBytes);
      offset += localHeader.length + nameBytes.length + entry.data.length;
    }

    const centralDirectory = concatBytes(centralParts);
    const endHeader = zipHeader(22, (view) => {
      view.setUint32(0, 0x06054b50, true);
      view.setUint16(4, 0, true);
      view.setUint16(6, 0, true);
      view.setUint16(8, entries.length, true);
      view.setUint16(10, entries.length, true);
      view.setUint32(12, centralDirectory.length, true);
      view.setUint32(16, offset, true);
      view.setUint16(20, 0, true);
    });
    return concatBytes([...fileParts, centralDirectory, endHeader]);
  }

  function downloadZip() {
    if (records.length === 0) return;

    const entries: ZipEntry[] = [
      { name: "nif_finder_results.tsv", data: textBytes(resultsTsv(records)) },
      { name: "nif_finder_results.csv", data: textBytes(resultsCsv(records)) },
    ];

    const nifFasta = fastaContent(records.filter((record) => nifGenes.includes(record.prediction)));
    if (nifFasta) {
      entries.push({ name: "nif_finder_detected_nif.faa", data: textBytes(nifFasta) });
    }

    const allHitFasta = fastaContent(records);
    if (allHitFasta) {
      entries.push({ name: "nif_finder_all_hits.faa", data: textBytes(allHitFasta) });
    }

    if (response?.plotPngBase64) {
      entries.push({ name: "nif_finder_scatter.png", data: base64Bytes(response.plotPngBase64) });
    }
    if (response?.genomicContextOverviewSvg) {
      entries.push({ name: "nif_finder_genome_overview.svg", data: textBytes(response.genomicContextOverviewSvg) });
    }
    if (response?.genomicContextLocalSvg) {
      entries.push({ name: "nif_finder_local_context.svg", data: textBytes(response.genomicContextLocalSvg) });
    }

    downloadBlob("nif_finder_results.zip", new Blob([createZip(entries)], { type: "application/zip" }));
  }

  function svgDataUri(svg: string) {
    return `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
  }

  async function analyze() {
    setLoading(true);
    setResponse(null);
    try {
      const requestBody = JSON.stringify({
        fasta,
        genbank: genbank.trim() || undefined,
        jobs,
        cpu,
        plot: plotOutput,
        evalue: Number(evalueThreshold),
      });
      let res: Response;
      if (requestBody.length > 3_500_000) {
        const tokenRes = await fetch("/api/compute-token", { method: "POST" });
        const tokenData = (await tokenRes.json()) as { apiUrl?: string; token?: string };
        if (tokenRes.ok && tokenData.apiUrl && tokenData.token) {
          res = await fetch(tokenData.apiUrl, {
            method: "POST",
            headers: {
              "content-type": "application/json",
              "x-analysis-token": tokenData.token,
            },
            body: requestBody,
          });
        } else {
          res = await fetch("/api/analyze", {
            method: "POST",
            headers: { "content-type": "application/json" },
            body: requestBody,
          });
        }
      } else {
        res = await fetch("/api/analyze", {
        method: "POST",
        headers: { "content-type": "application/json" },
          body: requestBody,
        });
      }
      const text = await res.text();
      const data = text
        ? (JSON.parse(text) as ApiResponse)
        : ({ error: `Analysis request returned an empty response (${res.status}).` } satisfies ApiResponse);
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

  function handleGenbankFile(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      setGenbank(String(reader.result ?? ""));
      setGenbankFileName(file.name);
    };
    reader.readAsText(file);
  }

  return (
    <main className="workspace">
      <aside className="sidebar">
        <div className="brand">
          <div>
            <h1>Nif-Finder</h1>
            <p>
              Web tool for detecting and classifying nitrogen fixation (<em>nif</em>) genes, including <em>nifH</em>,{" "}
              <em>nifD</em>, <em>nifK</em>, <em>nifE</em>, <em>nifN</em>, and <em>nifB</em>, from protein or genome FASTA
              using HMMER3 hmmscan and nearest-neighbour (1-NN) classification on homology and protein length plots.
            </p>
          </div>
        </div>

        <label className="field">
          Protein FASTA
          <textarea
            value={fasta}
            onChange={(event) => setFasta(event.target.value)}
            spellCheck={false}
          />
        </label>

        <label className="field compact-field example-field">
          Example dataset
          <select value={exampleDataset} onChange={(event) => void loadExampleDataset(event.target.value)}>
            {exampleDatasets.map((dataset) => (
              <option key={dataset.id} value={dataset.id}>
                {dataset.label}
              </option>
            ))}
          </select>
        </label>

        <div className="upload-row">
          <label className="file-button" title="Load protein FASTA from a local file">
            <FileUp size={18} aria-hidden />
            <span>Load protein FASTA</span>
            <input type="file" accept=".faa,.fa,.fasta,.txt" onChange={handleFile} />
          </label>
        </div>
        <p className="input-note">Protein FASTA max size: 10 MB.</p>

        <div className="upload-row genbank-upload">
          <label className="file-button" title="Load an optional GenBank file for genomic context plots">
            <FileUp size={18} aria-hidden />
            <span>Load GenBank</span>
            <input type="file" accept=".gb,.gbk,.gbff,.genbank,.txt" onChange={handleGenbankFile} />
          </label>
          {genbank ? (
            <button
              className="ghost-button compact-button"
              type="button"
              onClick={() => {
                setGenbank("");
                setGenbankFileName("");
              }}
            >
              Clear
            </button>
          ) : null}
        </div>
        <p className="input-note">
          Optional GenBank file for genome plotting. Max size: 25 MB.{" "}
          {genbankFileName ? `Loaded: ${genbankFileName}` : "No GenBank loaded."}
        </p>

        <div className="settings run-settings">
          <div className="section-title">
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
          <label>
            E-value threshold
            <input
              type="text"
              inputMode="decimal"
              value={evalueThreshold}
              onChange={(event) => setEvalueThreshold(event.target.value)}
            />
          </label>
        </div>
        <p className="input-note">E-value can affect sensitivity; the default is recommended.</p>

        <div className="settings advanced-settings">
          <div className="section-title">
            <Settings size={17} aria-hidden />
            Output
          </div>
          <label className="toggle-row">
            <input type="checkbox" checked={plotOutput} onChange={(event) => setPlotOutput(event.target.checked)} />
            Plot output
          </label>
          <label className="toggle-row">
            <input
              type="checkbox"
              checked={showOnlyNifHits}
              onChange={(event) => setShowOnlyNifHits(event.target.checked)}
            />
            Show only nif hits
          </label>
        </div>

        <button className="run-button" type="button" onClick={analyze} disabled={loading}>
          <Play size={18} aria-hidden />
          {loading ? "Running" : "Run analysis"}
        </button>
      </aside>

      <section className="results">
        <div className="summary-row">
          <div className="counter">{totalNifCopies} nif copies identified</div>
        </div>

        {response?.error ? (
          <div className="notice">
            <AlertCircle size={20} aria-hidden />
            <div>
              <strong>{response.error}</strong>
              {response.setup ? <p>{response.setup}</p> : null}
              {response.detail ? <p>{response.detail}</p> : null}
              {response.sequenceCount ? <p>{response.sequenceCount} FASTA record(s) were parsed.</p> : null}
            </div>
          </div>
        ) : null}

        {records.length > 0 ? (
          <>
            <div className="summary-table-wrap">
              <p className="table-caption">
                Table 1. Summary of <em>nif</em> genes identified by Nif-Finder.
              </p>
              <table className="summary-table">
                <thead>
                  <tr>
                    <th>Gene</th>
                    <th>Total copies</th>
                    <th>Full-length</th>
                    <th>Operon</th>
                    <th>Fragment</th>
                  </tr>
                </thead>
                <tbody>
                  {nifSummary.map((row) => (
                    <tr key={row.gene}>
                      <th scope="row">{row.gene}</th>
                      <td>{row.total}</td>
                      <td>{row.full}</td>
                      <td>{row.operon}</td>
                      <td>{row.incomplete}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            <div className="chart-panel">
              {response?.plotPngBase64 ? (
                <>
                  <img
                    className="nif-plot"
                    src={`data:image/png;base64,${response.plotPngBase64}`}
                    alt="Nif-Finder scatter plot"
                  />
                  <p className="figure-caption">{figure1Caption}</p>
                </>
              ) : (
                <div className="empty-state compact">No plot image was returned by the compute API.</div>
              )}
            </div>

            <div className="table-wrap">
              <p className="table-caption">{table2Caption}</p>
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
                  {displayedRecords.map((record) => (
                    <tr key={`${record.query}-${record.prediction}`}>
                      <td>{record.query}</td>
                      <td>{displayPrediction(record)}</td>
                      <td>{displayCompleteness(record)}</td>
                      <td>{record.logEvalue.toFixed(2)}</td>
                      <td>{record.alignLength}</td>
                      <td>{record.queryLength}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            {response?.genomicContextOverviewSvg || response?.genomicContextLocalSvg || response?.genomicContextMessage ? (
              <div className="genomic-context">
                {response?.genomicContextMessage ? (
                  <div className="notice">
                    <AlertCircle size={20} aria-hidden />
                    <div>
                      <strong>Genomic context was not plotted.</strong>
                      <p>{response.genomicContextMessage}</p>
                    </div>
                  </div>
                ) : null}

                {response?.genomicContextOverviewSvg ? (
                  <div className="chart-panel">
                    <img
                      className="context-plot"
                      src={svgDataUri(response.genomicContextOverviewSvg)}
                      alt="Whole genome view of matched nif gene locations"
                    />
                  </div>
                ) : null}

                {response?.genomicContextLocalSvg ? (
                  <div className="chart-panel">
                    <img
                      className="context-plot"
                      src={svgDataUri(response.genomicContextLocalSvg)}
                      alt="Enlarged genomic view around matched nif hits"
                    />
                  </div>
                ) : null}

                {response?.genomicContextOverviewSvg || response?.genomicContextLocalSvg ? (
                  <p className="figure-caption">{figure2Caption}</p>
                ) : null}
              </div>
            ) : null}

            <div className="download-row final-download-row" aria-label="Download results">
              <button className="ghost-button" type="button" onClick={downloadZip} disabled={records.length === 0}>
                <Download size={16} aria-hidden />
                Download ZIP
              </button>
            </div>
          </>
        ) : (
          <div className="empty-state">
            Submit a protein FASTA file to view predicted nif hits, completeness calls, and scatter plot output.
          </div>
        )}

        <footer className="site-footer">
          <p>
            Citation
            <br />
            Uesaka K, Fujita Y. Accurate prediction of nitrogen fixation in cyanobacteria reveals the dynamic evolution
            driving high retention rate with mosaic distribution. <em>bioRxiv</em>. 2026.{" "}
            <a href="https://doi.org/10.64898/2026.01.15.699626" target="_blank" rel="noreferrer">
              doi:10.64898/2026.01.15.699626
            </a>
          </p>
          <p>
            Code:{" "}
            <a href="https://github.com/kazumaxneo/Nif_finder" target="_blank" rel="noreferrer">
              github.com/kazumaxneo/Nif_finder
            </a>
          </p>
          <div className="software-citations">
            <p>Software citations:</p>
            <ul>
              <li>
                Eddy SR. Accelerated Profile HMM Searches. <em>PLoS Computational Biology</em>. 2011.{" "}
                <a href="https://doi.org/10.1371/journal.pcbi.1002195" target="_blank" rel="noreferrer">
                  doi:10.1371/journal.pcbi.1002195
                </a>
              </li>
              <li>
                Shimoyama Y. pyGenomeViz: A genome visualization python package for comparative genomics.{" "}
                <a href="https://github.com/moshi4/pyGenomeViz" target="_blank" rel="noreferrer">
                  github.com/moshi4/pyGenomeViz
                </a>
              </li>
            </ul>
          </div>
        </footer>
      </section>
    </main>
  );
}
