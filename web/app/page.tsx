"use client";

import { ChangeEvent, useMemo, useState } from "react";
import { Activity, AlertCircle, BarChart3, FileUp, Play, Settings } from "lucide-react";

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
  plotPngBase64?: string | null;
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

export default function Home() {
  const [fasta, setFasta] = useState("");
  const [jobs, setJobs] = useState(3);
  const [cpu, setCpu] = useState(6);
  const [loading, setLoading] = useState(false);
  const [response, setResponse] = useState<ApiResponse | null>(null);

  const records = response?.records ?? [];
  const nifSummary = useMemo(() => {
    return nifGenes.map((gene) => {
      const geneRecords = records.filter((record) => record.prediction === gene);
      const full = geneRecords.filter(
        (record) => record.completeness === "Full" || record.completeness === "Full_operon",
      ).length;
      const incomplete = geneRecords.filter((record) => record.completeness === "Fragment").length;

      return {
        gene,
        total: geneRecords.length,
        full,
        incomplete,
      };
    });
  }, [records]);
  const totalNifCopies = nifSummary.reduce((sum, row) => sum + row.total, 0);

  async function analyze() {
    setLoading(true);
    setResponse(null);
    try {
      const res = await fetch("/api/analyze", {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ fasta, jobs, cpu }),
      });
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
          <textarea
            value={fasta}
            onChange={(event) => setFasta(event.target.value)}
            placeholder="Paste protein FASTA here, or load a .faa/.fa/.fasta file."
            spellCheck={false}
          />
        </label>

        <div className="upload-row">
          <label className="file-button" title="Load protein FASTA from a local file">
            <FileUp size={18} aria-hidden />
            <span>Load protein FASTA</span>
            <input type="file" accept=".faa,.fa,.fasta,.txt" onChange={handleFile} />
          </label>
          <button className="ghost-button" type="button" onClick={() => setFasta(sampleFasta)}>
            Sample
          </button>
        </div>
        <p className="input-note">Max size: 10 MB. Sample is a demo set of known nifHDKENB proteins.</p>

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
          <div className="counter">{totalNifCopies} nif copies</div>
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
                    <th>Incomplete</th>
                  </tr>
                </thead>
                <tbody>
                  {nifSummary.map((row) => (
                    <tr key={row.gene}>
                      <th scope="row">{row.gene}</th>
                      <td>{row.total}</td>
                      <td>{row.full}</td>
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
