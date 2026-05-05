"use client";

import { ChangeEvent, useMemo, useState } from "react";
import { AlertCircle, BarChart3, Download, FileUp, Play, Settings } from "lucide-react";

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

const sampleFasta = `>Calothrix_sp.NIES-4101_01433|nifH
MSQDNIRQIAFYGKGGIGKSTTSQNTLAAMAEMGQRILIVGCDPKADSTRLMLHSKAQTT
VLHLAAERGAVEDIELEEVMLTGFRGVRCVESGGPEPGVGCAGRGIITAINFLEENGAYQ
DLDFVSYDVLGDVVCGGFAMPIREGKAQEIYIVTSGEMMAMYAANNIARGVLKYAHSGGV
RLGGLICNSRNTDREIELIETLASRLNTKMIHYVPRDNIVQHAELRRMTVNEYAPESNQS
NEYRTLAKKIINNKELDIPTPIEMEELEELLIEFGILESEENAAKLIASQDAAMKK
>Calothrix_sp.NIES-4101_03558|nifN
MAIVSSPNKPVVVNPLKQSQTLGAAMAFLGLKGMVPLLHGSQGCTAFAKVILVKHFREAI
PLATTAMTEVTTILGGEDNITQAILTLMEKSKPEIIGLCTTGLTETRGDDMQGILKDFRK
LHPELDSLPIILVSSPDFKGALQDGFAAAVEGIVRDVPEAGRTIPQQVNILAGSAFTPGD
IQQVKDIVSEFGLTPIVVPDLSASMDGHLDDSWSPITVGGTTLEQLKQVGSSAFTIALGE
SVRGAAEILEQRFEIPFEVFTELTGLQPVDEFLQALSDVSGVSVPQKFRHQRKQLQDAML
DTHFYFGRKRISLALEPDLLWSTVLFLQSMGAEIQAAVTTTRSPLLEQLPIDKVVIGDYE
DFEELAKGSDLLIGNSHAANIAKRLGIPLYRQGIPIIDRLGNGQFTKVGYKGTMDILFDI
GNILLETEEMKARDLKYEI
>Calothrix_sp.NIES-4101_03559|nifE
MKITQAKINELLTEPGCEHNSQKDSQKKKVCSQQAQPGAAQGGCAFDGAMIALVPITDAA
HLVHGPIACAGNSWGSRGSLSSGPTLYKMGFTTDLSENDVIFGGEKKLYKAILDIHKRYN
PAAIFVYSTCVTALIGEDMDAVCQAAAKKINTPVVPVNAAGFVGSKNLGNRIGGEALLEY
VVGTGEPDYTTPYDINLIGEYNIAGEMWNILPLLKKLGIRVLSKITGDARYQEITYAHRA
KLNVMICSKALINMARKMEERYGIPYIEESFYGVEDMNRCLRNIAAKLGDVELQARTEQL
IAEENAKLDVALAPYRARLKGRKVVIYTGGVKSWSIISAAQDLGMDVVATSTKKSTEEDK
AKIKKLLGKDGITLEKGNPQEILRVINETKAEMLIAGGRNQYTALKARIPFLDINQERHH
TYAGYSGMVEMARELDEALHSPVWEQIRKPAPWDFGTEILAESSLDWLDLPEEIFV
>Calothrix_sp.NIES-4101_03560|nifK
MFTEPVDFFIGNSYGKYLWRDTKIPMVRIGYPLFDRHHLHRYATLGYQGGLNLLNWVVNT
LLDEMDRNSNITGVTDISFDLIR
>Calothrix_sp.NIES-4101_03569|nifK
MPQNPDRIVDHLDLFKQEEYTKLFKNKRDNFEGAHSDAEVERVSEWTKSWEYREKNFARE
ALTVNPAKGCQPVGAMFAALGIEGTLPFVQGSQGCVAYFRTHLSRHYKEPCSAVSSSMTE
DAAVFGGLNNMIEGLEVSYKLYKPKMIAVCTTCMAEVIGDDLGAFITNAKGAGAIPQELP
VPFAHTPSFVGSHTTGYDNMMKGILTNLTEGKKGKANGKINVISGFDTYVGNNREVKRML
GLMGVEHTIVSDTSDYFDAPNDGEYQMYPGGTKIEDVADTINAKATVCLQAYTTQKTREY
IEKTWKQTTKVVRPFGVRGTDEFLTTISELTGKAIPAELEVERGRLVDAMTDSYAWIHGK
KFAIYGDPDLIISVTSFLLEMGAEPVHILCNNGDQEFKAEMEAMLKASPFGQGATVWIQK
DLWHFRRVICLLSN
>Calothrix_sp.NIES-4101_03641|nifD
MLLEEMGLRVVAQWSGDGTINELIQGPAAKLILIHCYRSMNYICRSLEESYGLPWMEFNF
FGPHKIAASLREIAAKFDKKIQDNAEKVIAKYQPTMDAVLKKYRPRLEGNTVMLYVGGLR
PRHVVPAFEDLGIKVVGTGYEFAHNDDYKRTTDYIDNATIIYDDVTAYEFEEFVKAKKPD
LIASGIKEKYVFQKMDYKYNVVRLISSMKYLKKPLCV
>Calothrix_sp.NIES-4101_03788|nifD
MTPPDKNLVEENQKLISEVLEAYPDKSKKKRQKHLNVHEEGKSDCGVKSNIKSVPGVMTA
RGCAYAGSKGVVWGPIKDMVHISHGPVGCGYWSWSGRRNYYVGKTGVDSFGTMHFTSDFQ
ERDIVFGGDKKLVKLIHEIDELFPLNKGISIQSECPIGLIGDDIEAVAKKTAKEIGKPVI
PLRCEGFRGVSQSLGHHIANDAIRDWVFPHFDKAKKENKLDFEPGPYDVALIGDYNISYI
GQHRKYRKVCLIL
>Calothrix_sp.NIES-4101_03795|nifB
MEEKIKERIAKHPCYSEEAHHHYARMHVAVAPACNIQCNYCNRKYDCANESRPGVVSELL
TPEEAAHKVLVVASKIPQMTVLGIAGPGDPLANPEKTFRTFELIADKAPDIKLCLSTNGL
MLPEYIDRIKQLNIDHVTITLNTIDPEIGAQIYSWVHYKRKRYRGVEGAKILLEKQMEGL
QALREADILCKVNSVMIPGINDQHLVEVNKMIRENGAFLHNIMPLISAPEHGTHFGLTGQ
RGPTSSELKSVQDSCSGNMKMMRHCRQCRADAVGLLGEDRSQEFTLSINCQNHPPPKFSA
SFGSQGLVYNPANGNR
>Calothrix_sp.NIES-4101_03823|nifH
MPIREGKAQEIYIVTSGEMMAMYAANNIARGILKYAHSGGVRLGGLICNSRKVDREAELI
ENLAERLNTQMIHFVPRDNIVQHAELRRMTVNEYAPDSNQSQEYRALAKKIINNEKLTIP
TPVDFTRGTYKTIERSSPPHNRVATPRRVSRLTYLHINKRAYIAPSLIFN
>Calothrix_sp.NIES-4101_03870|nifB
MDNCSDSDSLYYRIDSVYSCSQEFTKDKFMDMAPEYDLETRQEVQEGIEKFRAELKAAKE
QIALTNGTNSQGTAVKERTKILVAVATKGGGLVNQHFGHAKEFQVYEVDGSEVKYVAHRR
VDHFCQGGYGEKATLDNIINAISDCKAVLVSKIGESPKEKLNTAGIEVVESYDVIEKVAL
DYYQEFTAH
>Calothrix_sp.NIES-4101_03874|nifH
MSDERIRQIAFYGKGGIGKSTTSQNTLAAMAEMGQRIMIVGCDPKADSSCFY
>Calothrix_sp.NIES-4101_03890|nifH
MLHSKAQTTVLHLAAERGAVEDLELEEVMLAGFRGVKCVESGGPEPGVGCAGRGIITAIN
FLEENGAYQDLDFVSYDVLGDVVCGGFANYR
`;

const nifGenes = ["nifH", "nifD", "nifK", "nifE", "nifN", "nifB"];
const exampleDatasets = [
  { id: "none", label: "None" },
  { id: "known-nif", label: "Known nifHDKENB demo" },
];

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

  function loadExampleDataset(dataset: string) {
    setExampleDataset(dataset);
    setResponse(null);
    if (dataset === "known-nif") {
      setFasta(sampleFasta);
    } else {
      setFasta("");
    }
    setGenbank("");
    setGenbankFileName("");
  }

  function downloadText(filename: string, content: string, type = "text/plain;charset=utf-8") {
    const blob = new Blob([content], { type });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    link.click();
    URL.revokeObjectURL(url);
  }

  function downloadTsv() {
    const header = ["Query", "-log_Evalue", "Align_Len", "Query_Length", "Prediction", "Completeness", "Operon"];
    const rows = displayedRecords.map((record) => [
      record.query,
      record.logEvalue.toFixed(2),
      String(record.alignLength),
      record.queryLength == null ? "N/A" : String(record.queryLength),
      record.prediction,
      record.completeness,
      record.operonLabel ?? "",
    ]);
    downloadText("nif_finder_results.tsv", [header, ...rows].map((row) => row.join("\t")).join("\n") + "\n");
  }

  function displayPrediction(record: ResultRecord) {
    return record.completeness === "Full_operon" && record.operonLabel ? record.operonLabel : record.prediction;
  }

  function displayCompleteness(record: ResultRecord) {
    return record.completeness === "Full_operon" ? "Operon" : record.completeness;
  }

  function downloadPlot() {
    if (!response?.plotPngBase64) return;
    const link = document.createElement("a");
    link.href = `data:image/png;base64,${response.plotPngBase64}`;
    link.download = "nif_finder_scatter.png";
    link.click();
  }

  function downloadSvg(filename: string, svg: string | null | undefined) {
    if (!svg) return;
    downloadText(filename, svg, "image/svg+xml;charset=utf-8");
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
            <h1>Nif-Finder Web</h1>
            <p>NifHDKENB identification from protein fasta</p>
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

        <div className="upload-row">
          <label className="file-button" title="Load protein FASTA from a local file">
            <FileUp size={18} aria-hidden />
            <span>Load protein FASTA</span>
            <input type="file" accept=".faa,.fa,.fasta,.txt" onChange={handleFile} />
          </label>
        </div>

        <div className="upload-row genbank-upload">
          <label className="file-button" title="Load an optional GenBank file for genomic context plots">
            <FileUp size={18} aria-hidden />
            <span>Load GenBank (necessary only for nif plotting on the genome)</span>
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
          Optional GenBank file for genome plotting. {genbankFileName ? `Loaded: ${genbankFileName}` : "No GenBank loaded."}
        </p>

        <label className="field compact-field example-field">
          Example dataset
          <select value={exampleDataset} onChange={(event) => loadExampleDataset(event.target.value)}>
            {exampleDatasets.map((dataset) => (
              <option key={dataset.id} value={dataset.id}>
                {dataset.label}
              </option>
            ))}
          </select>
        </label>
        <p className="input-note">Max size: 10 MB. Example is a demo set of known nifHDKENB proteins.</p>

        <div className="settings run-settings">
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
        <p className="input-note">Changing E-value can affect sensitivity; the default is recommended.</p>

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
          <div>
            <p className="eyebrow">Results</p>
            <h2>Nif hit overview</h2>
          </div>
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
              <div className="section-title">
                <BarChart3 size={18} aria-hidden />
                Nif-Finder reference plot
              </div>
              {response?.plotPngBase64 ? (
                <img
                  className="nif-plot"
                  src={`data:image/png;base64,${response.plotPngBase64}`}
                  alt="Nif-Finder scatter plot"
                />
              ) : (
                <div className="empty-state compact">No plot image was returned by the compute API.</div>
              )}
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

            <div className="download-row" aria-label="Download results">
              <button className="ghost-button" type="button" onClick={downloadTsv} disabled={displayedRecords.length === 0}>
                <Download size={16} aria-hidden />
                Download TSV
              </button>
              <button className="ghost-button" type="button" onClick={downloadPlot} disabled={!response?.plotPngBase64}>
                <Download size={16} aria-hidden />
                Download PNG
              </button>
            </div>

            {response?.genomicContextOverviewSvg || response?.genomicContextLocalSvg || response?.genomicContextMessage ? (
              <div className="genomic-context">
                <div className="summary-row genomic-context-header">
                  <div>
                    <p className="eyebrow">GenBank context</p>
                    <h2>Nif gene locations</h2>
                  </div>
                  <div className="download-row" aria-label="Download genomic context">
                    <button
                      className="ghost-button"
                      type="button"
                      onClick={() => downloadSvg("nif_finder_genome_overview.svg", response?.genomicContextOverviewSvg)}
                      disabled={!response?.genomicContextOverviewSvg}
                    >
                      <Download size={16} aria-hidden />
                      Whole Genome SVG
                    </button>
                    <button
                      className="ghost-button"
                      type="button"
                      onClick={() => downloadSvg("nif_finder_local_context.svg", response?.genomicContextLocalSvg)}
                      disabled={!response?.genomicContextLocalSvg}
                    >
                      <Download size={16} aria-hidden />
                      Enlarged View SVG
                    </button>
                  </div>
                </div>

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
                    <div className="section-title">
                      <BarChart3 size={18} aria-hidden />
                      Whole Genome View
                    </div>
                    <img
                      className="context-plot"
                      src={svgDataUri(response.genomicContextOverviewSvg)}
                      alt="Whole genome view of matched nif gene locations"
                    />
                  </div>
                ) : null}

                {response?.genomicContextLocalSvg ? (
                  <div className="chart-panel">
                    <div className="section-title">
                      <BarChart3 size={18} aria-hidden />
                      Enlarged View
                    </div>
                    <img
                      className="context-plot"
                      src={svgDataUri(response.genomicContextLocalSvg)}
                      alt="Enlarged genomic view around matched nif hits"
                    />
                  </div>
                ) : null}
              </div>
            ) : null}
          </>
        ) : (
          <div className="empty-state">
            Submit a protein FASTA file to view predicted nif hits, completeness calls, and scatter plot output.
          </div>
        )}

        <footer className="site-footer">
          <p>
            Citation: Uesaka K, Fujita Y. Accurate prediction of nitrogen fixation in cyanobacteria reveals the
            dynamic evolution driving high retention rate with mosaic distribution. <em>bioRxiv</em>. 2026.{" "}
            <a href="https://doi.org/10.64898/2026.01.15.699626" target="_blank" rel="noreferrer">
              doi:10.64898/2026.01.15.699626
            </a>
          </p>
          <p>
            Source code:{" "}
            <a href="https://github.com/kazumaxneo/Nif_finder" target="_blank" rel="noreferrer">
              github.com/kazumaxneo/Nif_finder
            </a>
          </p>
        </footer>
      </section>
    </main>
  );
}
