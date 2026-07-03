"use client";

import { ChangeEvent, Dispatch, SetStateAction, useEffect, useMemo, useState } from "react";
import { AlertCircle, BookOpen, Download, FileText, FileUp, GitCompareArrows, HomeIcon, Info, Play } from "lucide-react";

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
  genomicContextGenbank?: string | null;
  genomicContextGenbankFilename?: string | null;
  vnfContextGenbank?: string | null;
  vnfContextGenbankFilename?: string | null;
  genomicContextMessage?: string | null;
  error?: string;
  setup?: string;
  detail?: string;
  sequenceCount?: number;
};

type VisitCountResponse = {
  enabled?: boolean;
  count?: number;
};

type ClusterRegion = {
  id: string;
  fileName: string;
  regionIndex: number;
  recordId: string;
  label: string;
  lengthBp: number;
  cdsCount: number;
  content: string;
};

type ClusterRegionResponse = {
  regions?: ClusterRegion[];
  warnings?: string[];
  error?: string;
  detail?: string;
};

type ClusterCompareResponse = {
  html?: string | null;
  plotFilename?: string | null;
  alignmentCsv?: string | null;
  warnings?: string[];
  error?: string;
  detail?: string;
};

type ClusterUploadFile = {
  id: string;
  slotIndex: number;
  name: string;
  uploadName: string;
  content: string;
};

type ClusterReference = {
  id: string;
  label: string;
};

let visitCountRequested = false;

const nifGenes = ["nifH", "nifD", "nifK", "nifE", "nifN", "nifB"];
const vnfGenes = ["vnfH/nifH", "vnfD", "vnfK", "vnfE/nifE", "vnfN/nifN", "vnfG", "vnfDG"];
const accessoryGenes = [
  "nifZ", "nifX", "nifP/cysE", "nifT", "nifV", "nifS", "nifU", "nifU_like",
  "modB/vupB", "modC/vupC", "modAlike", "vupA/modA", "vupB/modB",
  "vupC/modC", "cnfR/patB", "cnfR/patB_like",
];
const targetGenes = [...nifGenes, ...vnfGenes, ...accessoryGenes];
const coreModelGenes = ["nifH", "nifD", "nifK", "nifE", "nifN", "nifB"];
const additionalNifModelGenes = [
  { id: "nifZ", label: "nifZ" },
  { id: "nifX", label: "nifX" },
  { id: "nifP", label: "nifP/cysE" },
  { id: "nifT", label: "nifT" },
  { id: "nifV", label: "nifV" },
  { id: "nifS", label: "nifS" },
  { id: "nifU", label: "nifU/nifU_like" },
  { id: "modA", label: "modAlike" },
  { id: "modB", label: "modB/vupB" },
  { id: "modC", label: "modC/vupC" },
  { id: "cnfR-patB", label: "cnfR/patB" },
];
const vnfVupModelGenes = [
  { id: "vnfH", label: "vnfH" },
  { id: "vnfD", label: "vnfD" },
  { id: "vnfK", label: "vnfK" },
  { id: "vnfE", label: "vnfE" },
  { id: "vnfN", label: "vnfN" },
  { id: "vnfG", label: "vnfG/vnfDG" },
  { id: "vupA", label: "vupA/modA" },
  { id: "vupB", label: "vupB/modB" },
  { id: "vupC", label: "vupC/modC" },
];
const defaultAdditionalNifModelGenes = additionalNifModelGenes.map((gene) => gene.id);
const defaultVnfVupModelGenes = vnfVupModelGenes.map((gene) => gene.id);
const maxJobs = 4;
const maxCpu = 12;
const maxContextPaddingKb = 30;
const maxClusterUploads = 5;
const clusterReferences: Record<"groupI" | "groupII", ClusterReference[]> = {
  groupI: [
    { id: "groupI-dg5", label: "Leptolyngbya boryana dg5 nif-cluster" },
    { id: "groupI-ims101", label: "Trichodesmium erythraeum IMS101 nif-cluster" },
    { id: "groupI-atcc29413-I", label: "Trichormus variabilis ATCC 29413 nif-cluster I" },
    { id: "groupI-atcc29413-II", label: "Trichormus variabilis ATCC 29413 nif-cluster II" },
  ],
  groupII: [
    { id: "groupII-sodalinema-ab48", label: "Sodalinema sp. AB48 soda lake nif-cluster" },
    { id: "groupII-phormidium-ccy1219", label: "Phormidium sp. CCY1219 nif-cluster" },
  ],
};
const defaultClusterReferenceIds = {
  groupI: clusterReferences.groupI.map((reference) => reference.id),
  groupII: clusterReferences.groupII.map((reference) => reference.id),
};
const directComputeRequestThreshold = 1_000_000;
const apiRouteRetryMaxBytes = 4_000_000;
const exampleDatasets = [
  { id: "none", label: "None" },
  { id: "leptolyngbya-boryana-dg5", label: "Leptolyngbya boryana dg5" },
  { id: "anabaena-variabilis-atcc-29413", label: "Anabaena variabilis ATCC 29413" },
  { id: "calothrix-fragmented", label: "Calothrix strain's fragmented nif genes" },
];
const demoClusterFigures: Record<string, { src: string; title: string; caption: string }> = {
  "leptolyngbya-boryana-dg5": {
    src: "/examples/leptolyngbya_boryana_dg5_nif_cluster.png",
    title: "Reference nif-cluster for Leptolyngbya boryana dg5",
    caption:
      "Expected nif-cluster organization for the Leptolyngbya boryana dg5 demo dataset.",
  },
  "anabaena-variabilis-atcc-29413": {
    src: "/examples/anabaena_variabilis_atcc_29413_nif_cluster.png",
    title: "Reference nif-cluster for Anabaena variabilis ATCC 29413",
    caption:
      "Expected nif-cluster organization for the Anabaena variabilis ATCC 29413 demo dataset.",
  },
};

const figure1Caption = (
  <>
    Figure 1. 2D Similarity Plot of homology search for Nif proteins encoded by <em>nifHDKENB</em>, supported
    VnfD/VnfG/VnfK/VnfDG proteins, and related targets of
    various bacterial genomes. Background reference points show the relationship between –log10(E-value) from HMMER3
    searches using <em>nif</em> HMM profiles and hit protein length across proteomes of 586 cyanobacterial strains and
    other bacterial strains. Query sequences submitted to Nif-Finder are overlaid and highlighted on this 2D Similarity
    Plot, and are classified by their best matching target profile. The color indicates the single best hit
    protein from SWISS-PROT sequences. Circle, triangle, and star plots represent hits from complete genome assemblies,
    draft genome assemblies, and Nif fusion proteins, respectively.
  </>
);

const figure2Caption = (
  <>
    Figure 2. Genomic locations of <em>nif</em>/<em>vnf</em> genes identified by Nif-Finder. Colored arrows indicate{" "}
    target genes identified by Nif-Finder, and gray arrows indicate neighboring coding sequences. Enlarged regional views show
    local genomic neighborhoods around target hits. A whole-genome view is shown when it provides additional
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

type ActiveTab = "run" | "compare" | "figure2" | "manual" | "about";

const navigationTabs: Array<{
  id: ActiveTab;
  label: string;
  icon: typeof HomeIcon;
}> = [
  { id: "run", label: "Run", icon: HomeIcon },
  { id: "compare", label: "nif-cluster comparison", icon: GitCompareArrows },
  { id: "figure2", label: "Cyanobacterial figures", icon: FileText },
  { id: "manual", label: "Manual", icon: BookOpen },
  { id: "about", label: "About", icon: Info },
];

export default function Home() {
  const [activeTab, setActiveTab] = useState<ActiveTab>("run");
  const [fasta, setFasta] = useState("");
  const [genbank, setGenbank] = useState("");
  const [genbankFileName, setGenbankFileName] = useState("");
  const [jobs, setJobs] = useState(1);
  const [cpu, setCpu] = useState(4);
  const [contextPaddingKb, setContextPaddingKb] = useState(10);
  const [plotOutput, setPlotOutput] = useState(true);
  const [includeAdditionalNif, setIncludeAdditionalNif] = useState(false);
  const [includeVnfVup, setIncludeVnfVup] = useState(false);
  const [saveVnfRegionGbk, setSaveVnfRegionGbk] = useState(false);
  const [skipAccessoryOnlyLocalContext, setSkipAccessoryOnlyLocalContext] = useState(true);
  const [selectedAdditionalNifGenes, setSelectedAdditionalNifGenes] = useState<string[]>(defaultAdditionalNifModelGenes);
  const [selectedVnfVupGenes, setSelectedVnfVupGenes] = useState<string[]>(defaultVnfVupModelGenes);
  const [showOnlyNifHits, setShowOnlyNifHits] = useState(false);
  const [exampleDataset, setExampleDataset] = useState("none");
  const [submittedExampleDataset, setSubmittedExampleDataset] = useState("none");
  const [evalueThreshold, setEvalueThreshold] = useState("1e-10");
  const [loading, setLoading] = useState(false);
  const [response, setResponse] = useState<ApiResponse | null>(null);
  const [visitCount, setVisitCount] = useState<number | null>(null);
  const [clusterGroup, setClusterGroup] = useState<"groupI" | "groupII">("groupI");
  const [selectedClusterTemplateIds, setSelectedClusterTemplateIds] = useState<string[]>(defaultClusterReferenceIds.groupI);
  const [clusterFiles, setClusterFiles] = useState<ClusterUploadFile[]>([]);
  const [clusterRegions, setClusterRegions] = useState<ClusterRegion[]>([]);
  const [selectedClusterRegionIds, setSelectedClusterRegionIds] = useState<Record<string, string>>({});
  const [clusterSlotWarnings, setClusterSlotWarnings] = useState<Record<number, string[]>>({});
  const [clusterLoading, setClusterLoading] = useState(false);
  const [clusterResult, setClusterResult] = useState<ClusterCompareResponse | null>(null);

  const records = response?.records ?? [];
  const displayedRecords = showOnlyNifHits
    ? records.filter((record) => targetGenes.includes(record.prediction))
    : records;
  const targetSummary = useMemo(() => {
    return targetGenes.map((gene) => {
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
  const totalTargetCopies = targetSummary.reduce((sum, row) => sum + row.total, 0);

  useEffect(() => {
    if (visitCountRequested) return;
    visitCountRequested = true;
    let cancelled = false;

    fetch("/api/visit-count", { method: "POST", cache: "no-store" })
      .then((res) => (res.ok ? res.json() : null))
      .then((data: VisitCountResponse | null) => {
        if (!cancelled && data?.enabled && typeof data.count === "number") {
          setVisitCount(data.count);
        }
      })
      .catch(() => {
        if (!cancelled) {
          setVisitCount(null);
        }
      });

    return () => {
      cancelled = true;
    };
  }, []);
  const demoClusterFigure = records.length > 0 ? demoClusterFigures[submittedExampleDataset] : undefined;
  const clusterSlots = useMemo(() => Array.from({ length: maxClusterUploads }, (_, index) => index), []);
  const clusterWarnings = useMemo(() => Object.values(clusterSlotWarnings).flat(), [clusterSlotWarnings]);
  const clusterReferenceOptions = clusterReferences[clusterGroup];
  const selectedClusterRegions = useMemo(() => {
    return clusterFiles
      .map((file) => {
        const selectedId = selectedClusterRegionIds[file.uploadName];
        const regions = clusterRegions.filter((region) => region.fileName === file.uploadName);
        return regions.find((region) => region.id === selectedId) ?? null;
      })
      .filter((region): region is ClusterRegion => Boolean(region));
  }, [clusterFiles, clusterRegions, selectedClusterRegionIds]);
  const allLoadedClusterFilesSelected = selectedClusterRegions.length === clusterFiles.length;
  const canRunClusterComparison =
    allLoadedClusterFilesSelected && selectedClusterTemplateIds.length + selectedClusterRegions.length >= 2;

  function clampNumber(value: number, min: number, max: number) {
    if (!Number.isFinite(value)) return min;
    return Math.min(max, Math.max(min, Math.trunc(value)));
  }

  function toggleModelGene(
    gene: string,
    checked: boolean,
    setGenes: Dispatch<SetStateAction<string[]>>,
  ) {
    setGenes((current) => {
      if (checked) {
        return current.includes(gene) ? current : [...current, gene];
      }
      return current.filter((item) => item !== gene);
    });
  }

  function changeClusterGroup(group: "groupI" | "groupII") {
    setClusterGroup(group);
    setSelectedClusterTemplateIds(defaultClusterReferenceIds[group]);
    setClusterResult(null);
  }

  function toggleClusterTemplate(templateId: string) {
    setSelectedClusterTemplateIds((current) =>
      current.includes(templateId) ? current.filter((item) => item !== templateId) : [...current, templateId],
    );
    setClusterResult(null);
  }

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
    } else if (dataset === "anabaena-variabilis-atcc-29413") {
      const [fastaResponse, genbankResponse] = await Promise.all([
        fetch("/examples/anabaena_variabilis_atcc_29413.faa"),
        fetch("/examples/anabaena_variabilis_atcc_29413.gbff"),
      ]);
      setFasta(await fastaResponse.text());
      setGenbank(await genbankResponse.text());
      setGenbankFileName("anabaena_variabilis_atcc_29413.gbff");
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

  function stableFileId(file: File, index: number) {
    return `${file.name}-${file.size}-${file.lastModified}-${index}`;
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

    const nifFasta = fastaContent(records.filter((record) => targetGenes.includes(record.prediction)));
    if (nifFasta) {
      entries.push({ name: "nif_finder_detected_nif_vnf_related.faa", data: textBytes(nifFasta) });
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
    const localContextGenbankFilename = response?.genomicContextGenbankFilename || "nif_finder_local_context.gbk";
    if (response?.genomicContextGenbank) {
      entries.push({ name: localContextGenbankFilename, data: textBytes(response.genomicContextGenbank) });
    }
    const vnfContextGenbankFilename = response?.vnfContextGenbankFilename || "nif_finder_vnf_context.gbk";
    if (response?.vnfContextGenbank) {
      entries.push({ name: vnfContextGenbankFilename, data: textBytes(response.vnfContextGenbank) });
    }

    downloadBlob("nif_finder_results.zip", new Blob([createZip(entries)], { type: "application/zip" }));
  }

  function downloadLocalContextGenbank() {
    if (!response?.genomicContextGenbank) return;
    downloadText(
      response.genomicContextGenbankFilename || "nif_finder_local_context.gbk",
      response.genomicContextGenbank,
      "chemical/x-genbank;charset=utf-8",
    );
  }

  function downloadVnfContextGenbank() {
    if (!response?.vnfContextGenbank) return;
    downloadText(
      response.vnfContextGenbankFilename || "nif_finder_vnf_context.gbk",
      response.vnfContextGenbank,
      "chemical/x-genbank;charset=utf-8",
    );
  }

  function svgDataUri(svg: string) {
    return `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
  }

  function summarizeNonJsonResponse(text: string) {
    return text.replace(/<[^>]*>/g, " ").replace(/\s+/g, " ").trim().slice(0, 600) || text.slice(0, 600);
  }

  async function readApiResponse(res: Response): Promise<{ data: ApiResponse; nonJson: boolean }> {
    const text = await res.text();
    if (!text) {
      return {
        data: { error: `Analysis request returned an empty response (${res.status}).` },
        nonJson: false,
      };
    }
    try {
      return { data: JSON.parse(text) as ApiResponse, nonJson: false };
    } catch {
      return {
        data: {
          error: `Analysis service returned a non-JSON response (${res.status}).`,
          detail: summarizeNonJsonResponse(text),
        },
        nonJson: true,
      };
    }
  }

  async function readJsonResponse<T extends { error?: string; detail?: string }>(res: Response): Promise<T> {
    const text = await res.text();
    if (!text) {
      return { error: `Request returned an empty response (${res.status}).` } as T;
    }
    try {
      return JSON.parse(text) as T;
    } catch {
      return {
        error: `Service returned a non-JSON response (${res.status}).`,
        detail: summarizeNonJsonResponse(text),
      } as T;
    }
  }

  async function postViaApiRoute(requestBody: string) {
    return fetch("/api/analyze", {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: requestBody,
    });
  }

  async function analyze() {
    setLoading(true);
    setResponse(null);
    setSubmittedExampleDataset(exampleDataset);
    try {
      const requestBody = JSON.stringify({
        fasta,
        genbank: genbank.trim() || undefined,
        jobs,
        cpu,
        contextPaddingKb,
        plot: plotOutput,
        evalue: Number(evalueThreshold),
        vnfMode: false,
        saveVnfRegionGbk: includeVnfVup && saveVnfRegionGbk,
        skipAccessoryOnlyLocalContext,
        selectedModelGenes: [
          ...(includeAdditionalNif ? selectedAdditionalNifGenes : []),
          ...(includeVnfVup ? selectedVnfVupGenes : []),
        ],
      });
      let res: Response;
      let usedDirectCompute = false;
      if (requestBody.length > directComputeRequestThreshold) {
        const tokenRes = await fetch("/api/compute-token", { method: "POST" });
        const { data: tokenData } = await readApiResponse(tokenRes) as { data: ApiResponse & { apiUrl?: string; token?: string }; nonJson: boolean };
        if (tokenRes.ok && tokenData.apiUrl && tokenData.token) {
          usedDirectCompute = true;
          res = await fetch(tokenData.apiUrl, {
            method: "POST",
            headers: {
              "content-type": "application/json",
              "x-analysis-token": tokenData.token,
            },
            body: requestBody,
          });
        } else {
          res = await postViaApiRoute(requestBody);
        }
      } else {
        res = await postViaApiRoute(requestBody);
      }
      let { data, nonJson } = await readApiResponse(res);
      if (usedDirectCompute && nonJson && requestBody.length <= apiRouteRetryMaxBytes) {
        res = await postViaApiRoute(requestBody);
        ({ data, nonJson } = await readApiResponse(res));
      }
      if (!res.ok && !data.error) {
        setResponse({
          ...data,
          error: data.detail || `Analysis request failed (${res.status}).`,
        });
      } else {
        setResponse(data);
      }
    } catch (error) {
      setResponse({ error: error instanceof Error ? error.message : "Request failed." });
    } finally {
      setLoading(false);
    }
  }

  function clusterUploadName(slotIndex: number, fileName: string) {
    return `slot_${slotIndex + 1}_${fileName.replace(/[^\w.-]+/g, "_")}`;
  }

  async function loadClusterFile(event: ChangeEvent<HTMLInputElement>, slotIndex: number) {
    const file = event.target.files?.[0];
    if (!file) return;
    setClusterLoading(true);
    setClusterResult(null);
    try {
      const uploadName = clusterUploadName(slotIndex, file.name);
      const loaded: ClusterUploadFile = {
        id: stableFileId(file, slotIndex),
        slotIndex,
        name: file.name,
        uploadName,
        content: await file.text(),
      };
      const res = await fetch("/api/cluster-regions", {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({
          files: [{ name: loaded.uploadName, content: loaded.content }],
        }),
      });
      const data = await readJsonResponse<ClusterRegionResponse>(res);
      if (!res.ok || data.error) {
        setClusterFiles((current) => current.filter((item) => item.slotIndex !== slotIndex));
        setClusterRegions((current) => current.filter((region) => region.fileName !== uploadName));
        setSelectedClusterRegionIds((current) => {
          const next = { ...current };
          delete next[uploadName];
          return next;
        });
        setClusterSlotWarnings((current) => ({
          ...current,
          [slotIndex]: [data.error ?? "Could not inspect GenBank regions.", data.detail ?? ""].filter(Boolean),
        }));
        return;
      }
      const regions = data.regions ?? [];
      setClusterFiles((current) =>
        [...current.filter((item) => item.slotIndex !== slotIndex), loaded].sort((a, b) => a.slotIndex - b.slotIndex),
      );
      setClusterRegions((current) => [
        ...current.filter((region) => region.fileName !== uploadName),
        ...regions,
      ]);
      setSelectedClusterRegionIds((current) => {
        const next = { ...current };
        delete next[uploadName];
        if (regions.length === 1) {
          next[uploadName] = regions[0].id;
        }
        return next;
      });
      setClusterSlotWarnings((current) => ({ ...current, [slotIndex]: data.warnings ?? [] }));
    } catch (error) {
      setClusterFiles((current) => current.filter((item) => item.slotIndex !== slotIndex));
      setClusterSlotWarnings((current) => ({
        ...current,
        [slotIndex]: [error instanceof Error ? error.message : "Could not read GenBank file."],
      }));
    } finally {
      setClusterLoading(false);
      event.target.value = "";
    }
  }

  function clearClusterFile(slotIndex: number) {
    const file = clusterFiles.find((item) => item.slotIndex === slotIndex);
    setClusterFiles((current) => current.filter((item) => item.slotIndex !== slotIndex));
    if (file) {
      setClusterRegions((current) => current.filter((region) => region.fileName !== file.uploadName));
      setSelectedClusterRegionIds((current) => {
        const next = { ...current };
        delete next[file.uploadName];
        return next;
      });
    }
    setClusterSlotWarnings((current) => {
      const next = { ...current };
      delete next[slotIndex];
      return next;
    });
    setClusterResult(null);
  }

  function clearClusterFiles() {
    setClusterFiles([]);
    setClusterRegions([]);
    setSelectedClusterRegionIds({});
    setClusterSlotWarnings({});
    setClusterResult(null);
  }

  async function runClusterComparison() {
    setClusterLoading(true);
    setClusterResult(null);
    try {
      const res = await fetch("/api/compare-clusters", {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({
          group: clusterGroup,
          templateIds: selectedClusterTemplateIds,
          regions: selectedClusterRegions,
        }),
      });
      const data = await readJsonResponse<ClusterCompareResponse>(res);
      setClusterResult(data);
    } catch (error) {
      setClusterResult({ error: error instanceof Error ? error.message : "Cluster comparison failed." });
    } finally {
      setClusterLoading(false);
    }
  }

  function downloadClusterHtml() {
    if (!clusterResult?.html) return;
    downloadText(clusterResult.plotFilename || "nif_cluster_comparison.html", clusterResult.html, "text/html;charset=utf-8");
  }

  function downloadClusterAlignment() {
    if (!clusterResult?.alignmentCsv) return;
    downloadText("nif_cluster_comparison_alignments.csv", clusterResult.alignmentCsv, "text/csv;charset=utf-8");
  }

  function handleFile(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    if (!file) return;
    setExampleDataset("none");
    setSubmittedExampleDataset("none");
    const reader = new FileReader();
    reader.onload = () => setFasta(String(reader.result ?? ""));
    reader.readAsText(file);
  }

  function handleGenbankFile(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    if (!file) return;
    setExampleDataset("none");
    setSubmittedExampleDataset("none");
    const reader = new FileReader();
    reader.onload = () => {
      setGenbank(String(reader.result ?? ""));
      setGenbankFileName(file.name);
    };
    reader.readAsText(file);
  }

  return (
    <main className="workspace">
      <header className="top-menu">
        <div className="top-menu-brand-row" aria-label="Nif-Finder site header">
          <span className="top-menu-name">Nif-Finder</span>
          <span className="top-menu-title">Nitrogen fixation gene finder</span>
        </div>
        <nav className="top-tabs" aria-label="Nif-Finder sections">
          {navigationTabs.map(({ id, label, icon: Icon }) => (
            <button
              key={id}
              className={activeTab === id ? "top-tab active" : "top-tab"}
              type="button"
              onClick={() => setActiveTab(id)}
            >
              <Icon aria-hidden="true" size={15} strokeWidth={2} />
              <span>{label}</span>
            </button>
          ))}
        </nav>
      </header>

      <aside className="sidebar">
        <div className="brand">
          <div className="brand-copy">
            <p>
              Web tool for detecting and classifying nitrogen fixation (<em>nif</em>/<em>vnf</em>) genes and related targets, including{" "}
              <em>nifH</em>, <em>nifD</em>, <em>nifK</em>, <em>nifE</em>, <em>nifN</em>, <em>nifB</em>,{" "}
              <em>vnfH/nifH</em>, <em>vnfD</em>, <em>vnfK</em>, <em>vnfE/nifE</em>, <em>vnfN/nifN</em>, <em>vnfG</em>, <em>vnfDG</em>, <em>nifZ</em>, <em>nifX</em>, <em>nifP/cysE</em>, <em>nifT</em>, <em>nifV</em>, <em>nifS</em>,{" "}
              <em>nifU</em>, <em>nifU_like</em>, <em>modAlike</em>, <em>modB/vupB</em>, <em>modC/vupC</em>, <em>vupA/modA</em>, <em>vupB/modB</em>, <em>vupC/modC</em>, and <em>cnfR/patB</em>, from protein or genome FASTA
              using HMMER3 hmmscan and nearest-neighbour (1-NN) classification on homology and protein length plots.
            </p>
            <div className="brand-citation">
              <p>
                Citation
                <br />
                Uesaka K, Fujita Y. Accurate prediction of nitrogen fixation in cyanobacteria reveals the dynamic
                evolution driving high retention rate with mosaic distribution. <em>bioRxiv</em>. 2026.{" "}
                <a href="https://doi.org/10.64898/2026.01.15.699626" target="_blank" rel="noreferrer">
                  doi:10.64898/2026.01.15.699626
                </a>
              </p>
              {activeTab === "run" ? (
                <>
                  <p className="brand-sub-site">
                    Sub site:{" "}
                    <a href="https://web-theta-black-17.vercel.app" target="_blank" rel="noreferrer">
                      https://web-theta-black-17.vercel.app
                    </a>
                  </p>
                  {visitCount !== null ? (
                    <p className="brand-visit-counter">Total visits: {visitCount.toLocaleString()}</p>
                  ) : null}
                </>
              ) : null}
            </div>
          </div>
          <img className="brand-figure" src="/nif_scatter_panels_transparent.png" alt="" aria-hidden="true" />
        </div>

        {activeTab === "run" ? (
          <>
        <label className="field">
          Protein FASTA
          <textarea
            value={fasta}
            onChange={(event) => {
              setFasta(event.target.value);
              setExampleDataset("none");
              setSubmittedExampleDataset("none");
            }}
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
          Optional GenBank file for genome plotting. Max size: 30 MB.{" "}
          {genbankFileName ? `Loaded: ${genbankFileName}` : "No GenBank loaded."}
        </p>

        <div className="settings-row">
          <div className="settings run-settings">
            <div className="section-title">
              Run settings
            </div>
            <label>
              Jobs
              <input
                type="number"
                min={1}
                max={maxJobs}
                value={jobs}
                onChange={(event) => setJobs(clampNumber(Number(event.target.value), 1, maxJobs))}
              />
            </label>
            <label>
              CPU
              <input
                type="number"
                min={1}
                max={maxCpu}
                value={cpu}
                onChange={(event) => setCpu(clampNumber(Number(event.target.value), 1, maxCpu))}
              />
            </label>
            <label className="evalue-field">
              E-value threshold
              <input
                type="text"
                inputMode="decimal"
                value={evalueThreshold}
                onChange={(event) => setEvalueThreshold(event.target.value)}
              />
            </label>
            <label className="context-padding-field">
              Size of the nif/vnf-related region to be visualized (kb)
              <input
                type="number"
                min={1}
                max={maxContextPaddingKb}
                value={contextPaddingKb}
                onChange={(event) => setContextPaddingKb(clampNumber(Number(event.target.value), 1, maxContextPaddingKb))}
              />
            </label>
          </div>

          <div className="settings advanced-settings">
            <div className="section-title">
              Output
            </div>
            <label className="toggle-row">
              <input type="checkbox" checked={plotOutput} onChange={(event) => setPlotOutput(event.target.checked)} />
              Plot output
            </label>
            <label className="toggle-row">
              <input
                type="checkbox"
                checked={includeAdditionalNif}
                onChange={(event) => setIncludeAdditionalNif(event.target.checked)}
              />
              Add additional nif-related detection
              <span className="experimental-label">(experimental)</span>
            </label>
            <label className="toggle-row">
              <input
                type="checkbox"
                checked={includeVnfVup}
                onChange={(event) => {
                  setIncludeVnfVup(event.target.checked);
                  if (!event.target.checked) {
                    setSaveVnfRegionGbk(false);
                  }
                }}
              />
              Add vnf/vup detection
              <span className="experimental-label">(experimental)</span>
            </label>
            <label className="toggle-row">
              <input
                type="checkbox"
                checked={showOnlyNifHits}
                onChange={(event) => setShowOnlyNifHits(event.target.checked)}
              />
              Show only target hits
            </label>
            <label className="toggle-row">
              <input
                type="checkbox"
                checked={skipAccessoryOnlyLocalContext}
                onChange={(event) => setSkipAccessoryOnlyLocalContext(event.target.checked)}
              />
              Hide local maps with only nifP/mod/vup hits
            </label>
          </div>
        </div>
        {includeAdditionalNif ? (
        <details className="experimental-options" open>
          <summary>Additional nif-related target genes <span>(experimental)</span></summary>
          <div className="experimental-options-body">
            <p>
              Core <em>nifHDKENB</em> models always run. Select additional available models to add to this analysis.
            </p>
            <div className="model-select-actions">
              <button type="button" className="ghost-button compact-button" onClick={() => setSelectedAdditionalNifGenes(defaultAdditionalNifModelGenes)}>
                All select
              </button>
              <button type="button" className="ghost-button compact-button" onClick={() => setSelectedAdditionalNifGenes([])}>
                Clear
              </button>
            </div>
            <div className="core-gene-list" aria-label="Core models">
              {coreModelGenes.map((gene) => (
                <span key={gene}>{gene}</span>
              ))}
            </div>
            <div className="model-checkbox-grid">
              {additionalNifModelGenes.map((gene) => (
                <label key={gene.id} className="toggle-row">
                  <input
                    type="checkbox"
                    checked={selectedAdditionalNifGenes.includes(gene.id)}
                    onChange={(event) => toggleModelGene(gene.id, event.target.checked, setSelectedAdditionalNifGenes)}
                  />
                  {gene.label}
                </label>
              ))}
            </div>
          </div>
        </details>
        ) : null}
        {includeVnfVup ? (
        <details className="experimental-options" open>
          <summary>Vnf/vup target genes <span>(experimental)</span></summary>
          <div className="experimental-options-body">
            <p>
              Core <em>nifHDKENB</em> models always run. Select vnf/vup models to add to this analysis.
            </p>
            <div className="model-select-actions">
              <button type="button" className="ghost-button compact-button" onClick={() => setSelectedVnfVupGenes(defaultVnfVupModelGenes)}>
                All select
              </button>
              <button type="button" className="ghost-button compact-button" onClick={() => setSelectedVnfVupGenes([])}>
                Clear
              </button>
            </div>
            <div className="core-gene-list" aria-label="Core models">
              {coreModelGenes.map((gene) => (
                <span key={gene}>{gene}</span>
              ))}
            </div>
            <div className="model-checkbox-grid">
              {vnfVupModelGenes.map((gene) => (
                <label key={gene.id} className="toggle-row">
                  <input
                    type="checkbox"
                    checked={selectedVnfVupGenes.includes(gene.id)}
                    onChange={(event) => toggleModelGene(gene.id, event.target.checked, setSelectedVnfVupGenes)}
                  />
                  {gene.label}
                </label>
              ))}
            </div>
            <label className="toggle-row">
              <input
                type="checkbox"
                checked={saveVnfRegionGbk}
                onChange={(event) => setSaveVnfRegionGbk(event.target.checked)}
              />
              Output vnf-only GenBank region
              <span className="experimental-label">(experimental)</span>
            </label>
          </div>
        </details>
        ) : null}
        <p className="input-note">
          E-value can affect sensitivity; the default is recommended.
        </p>

        <button className="run-button" type="button" onClick={analyze} disabled={loading}>
          <Play size={18} aria-hidden />
          {loading ? "Running" : "Run analysis"}
        </button>
        {loading ? (
          <p className="input-note running-note">Analysis is running. Results may take several minutes</p>
        ) : null}
          </>
        ) : null}

        {activeTab === "compare" ? (
          <div className="cluster-controls">
            <div className="section-title">nif-cluster comparison</div>
            <label className="field compact-field cluster-group-field">
              Comparison group
              <select value={clusterGroup} onChange={(event) => changeClusterGroup(event.target.value as "groupI" | "groupII")}>
                <option value="groupI">Group I nif</option>
                <option value="groupII">Group II nif</option>
              </select>
            </label>
            <p className="input-note cluster-note">
              Select reference clusters to include, then optionally upload up to 5 Nif-Finder GenBank regions.
              At least two total regions are required for comparison; each selected region must be 100 kb or smaller.
            </p>
            <div className="cluster-reference-options" aria-label="Reference clusters">
              {clusterReferenceOptions.map((reference) => (
                <label className="cluster-reference-option" key={reference.id}>
                  <input
                    type="checkbox"
                    checked={selectedClusterTemplateIds.includes(reference.id)}
                    onChange={() => toggleClusterTemplate(reference.id)}
                  />
                  <span>{reference.label}</span>
                </label>
              ))}
            </div>
            <div className="cluster-slot-list">
              {clusterSlots.map((slotIndex) => {
                const file = clusterFiles.find((item) => item.slotIndex === slotIndex);
                const regions = file ? clusterRegions.filter((region) => region.fileName === file.uploadName) : [];
                return (
                  <div className="cluster-slot-row" key={slotIndex}>
                    <span className="cluster-slot-number">{slotIndex + 1}</span>
                    <label className="file-button cluster-file-button" title="Load a Nif-Finder GenBank region for clinker comparison">
                      <FileUp size={18} aria-hidden />
                      <span>GBK</span>
                      <input type="file" accept=".gb,.gbk,.gbff,.genbank,.txt" onChange={(event) => loadClusterFile(event, slotIndex)} />
                    </label>
                    <span className="cluster-file-name" title={file?.name ?? ""}>
                      {file ? file.name : "No file"}
                    </span>
                    <select
                      className="cluster-region-select"
                      value={file ? selectedClusterRegionIds[file.uploadName] ?? "" : ""}
                      onChange={(event) =>
                        file
                          ? setSelectedClusterRegionIds((current) => ({ ...current, [file.uploadName]: event.target.value }))
                          : undefined
                      }
                      disabled={!file || regions.length === 0}
                      aria-label={`Select region for cluster GBK ${slotIndex + 1}`}
                    >
                      <option value="">
                        {file ? "Select one region" : "No file loaded"}
                      </option>
                      {regions.map((region) => (
                        <option key={region.id} value={region.id}>
                          region {region.regionIndex}: {region.recordId} ({region.lengthBp.toLocaleString()} bp, {region.cdsCount} CDS)
                        </option>
                      ))}
                    </select>
                    <button
                      className="ghost-button compact-button cluster-clear-button"
                      type="button"
                      onClick={() => clearClusterFile(slotIndex)}
                      disabled={!file}
                    >
                      Clear
                    </button>
                  </div>
                );
              })}
            </div>
            {clusterFiles.length > 0 ? (
              <button className="ghost-button compact-button cluster-clear-all-button" type="button" onClick={clearClusterFiles}>
                Clear all
              </button>
            ) : null}
            <p className="input-note">
              {clusterFiles.length > 0 ? `${clusterFiles.length} file(s) loaded.` : "No cluster GenBank loaded."}
            </p>
            {clusterWarnings.length > 0 ? (
              <div className="notice compact-notice">
                <AlertCircle size={18} aria-hidden />
                <div>
                  <strong>Cluster input notice.</strong>
                  {clusterWarnings.map((warning) => (
                    <p key={warning}>{warning}</p>
                  ))}
                </div>
              </div>
            ) : null}
            <button
              className="run-button"
              type="button"
              onClick={runClusterComparison}
              disabled={clusterLoading || !canRunClusterComparison}
            >
              <GitCompareArrows size={18} aria-hidden />
              {clusterLoading ? "Running" : "Run clinker"}
            </button>
          </div>
        ) : null}
      </aside>

      {activeTab === "run" ? (
      <section className="results">
        <div className="summary-row">
          <div className="counter">{totalTargetCopies} target copies identified</div>
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
                Table 1. Summary of target genes identified by Nif-Finder.
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
                  {targetSummary.map((row) => (
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

            {response?.genomicContextOverviewSvg || response?.genomicContextLocalSvg || response?.genomicContextMessage ? (
              <div className="genomic-context">
                {response?.genomicContextMessage ? (
                  <div className="notice">
                    <AlertCircle size={20} aria-hidden />
                    <div>
                      <strong>Genomic context notice.</strong>
                      <p>{response.genomicContextMessage}</p>
                    </div>
                  </div>
                ) : null}

                {response?.genomicContextOverviewSvg ? (
                  <div className="chart-panel">
                    <img
                      className="context-plot"
                      src={svgDataUri(response.genomicContextOverviewSvg)}
                      alt="Whole genome view of matched target gene locations"
                    />
                  </div>
                ) : null}

                {response?.genomicContextLocalSvg ? (
                  <div className="chart-panel">
                    <img
                      className="context-plot"
                      src={svgDataUri(response.genomicContextLocalSvg)}
                      alt="Enlarged genomic view around matched target hits"
                    />
                  </div>
                ) : null}

                {response?.genomicContextOverviewSvg || response?.genomicContextLocalSvg ? (
                  <p className="figure-caption">{figure2Caption}</p>
                ) : null}

                {response?.genomicContextGenbank ? (
                  <div className="download-row" aria-label="Download Figure 2 genomic context">
                    <button className="ghost-button" type="button" onClick={downloadLocalContextGenbank}>
                      <Download size={16} aria-hidden />
                      Download nif/vnf-related region (gbk)
                    </button>
                  </div>
                ) : null}
              </div>
            ) : null}

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

            {demoClusterFigure ? (
              <div className="chart-panel demo-cluster-panel">
                <h3>{demoClusterFigure.title}</h3>
                <img
                  className="demo-cluster-plot"
                  src={demoClusterFigure.src}
                  alt={demoClusterFigure.title}
                />
                <p className="figure-caption">{demoClusterFigure.caption}</p>
              </div>
            ) : null}

            <div className="download-row final-download-row" aria-label="Download results">
              {response?.vnfContextGenbank ? (
                <button className="ghost-button" type="button" onClick={downloadVnfContextGenbank}>
                  <Download size={16} aria-hidden />
                  Download vnf region (gbk)
                </button>
              ) : null}
              <button className="ghost-button" type="button" onClick={downloadZip} disabled={records.length === 0}>
                <Download size={16} aria-hidden />
                Download ZIP
              </button>
            </div>
          </>
        ) : (
          <div className="empty-state">
            Submit a protein FASTA file to view predicted target hits, completeness calls, and scatter plot output.
          </div>
        )}

        <footer className="site-footer">
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
      ) : (
        <section className={activeTab === "compare" && clusterResult?.html ? "info-page wide-clinker-page" : "info-page"}>
          {activeTab === "compare" ? (
            <article className="manual-body cluster-comparison-page">
              <h2>nif-cluster comparison</h2>
              <figure className="cluster-reference-figure">
                <figcaption>Typical gene organization of Group I <em>nif</em>-clusters.</figcaption>
                <img src="/group-i-nif-cluster.jpg" alt="Group I nif-cluster reference gene maps" />
              </figure>
              <p>
                Compare Nif-Finder GenBank cluster regions with clinker. For Group I nif, the dg5 and ATCC29413
                teaching clusters are added automatically; uploaded regions are compared against those references.
              </p>

              {clusterResult?.error ? (
                <div className="notice">
                  <AlertCircle size={20} aria-hidden />
                  <div>
                    <strong>{clusterResult.error}</strong>
                    {clusterResult.detail ? <p>{clusterResult.detail}</p> : null}
                  </div>
                </div>
              ) : null}

              {clusterResult?.warnings?.length ? (
                <div className="notice">
                  <AlertCircle size={20} aria-hidden />
                  <div>
                    <strong>Comparison notice.</strong>
                    {clusterResult.warnings.map((warning) => (
                      <p key={warning}>{warning}</p>
                    ))}
                  </div>
                </div>
              ) : null}

              {clusterResult?.html ? (
                <>
                  <div className="download-row">
                    <button className="ghost-button" type="button" onClick={downloadClusterHtml}>
                      <Download size={16} aria-hidden />
                      Download interactive clinker HTML
                    </button>
                    {clusterResult.alignmentCsv ? (
                      <button className="ghost-button" type="button" onClick={downloadClusterAlignment}>
                        <Download size={16} aria-hidden />
                        Download alignments CSV
                      </button>
                    ) : null}
                  </div>
                  <p className="input-note cluster-download-note">
                    The embedded clinker Save SVG button is enabled. If your browser still blocks it, download the
                    interactive clinker HTML and open it directly, then use Save SVG there.
                  </p>
                  <iframe
                    className="clinker-frame"
                    title="clinker nif-cluster comparison"
                    sandbox="allow-scripts allow-downloads"
                    srcDoc={clusterResult.html}
                  />
                  <div className="software-citations cluster-citation">
                    <p>clinker citation:</p>
                    <ul>
                      <li>
                        Gilchrist CLM, Chooi YH. clinker &amp; clustermap.js: automatic generation of gene cluster
                        comparison figures. <em>Bioinformatics</em>. 2021.{" "}
                        <a href="https://doi.org/10.1093/bioinformatics/btab007" target="_blank" rel="noreferrer">
                          doi:10.1093/bioinformatics/btab007
                        </a>
                      </li>
                    </ul>
                  </div>
                </>
              ) : (
                <div className="empty-state">
                  Upload up to 5 Nif-Finder GenBank region files, select one region per file when needed, and run clinker.
                </div>
              )}
            </article>
          ) : null}

          {activeTab === "figure2" ? (
            <article className="manual-body figure-page">
              <div className="revision-note">
                <p>Kazuma Uesaka and Yuichi Fujita, 2026.</p>
                <p>
                  &quot;Accurate prediction of nitrogen fixation in cyanobacteria reveals the dynamic evolution driving
                  high retention rate with mosaic distribution&quot;
                </p>
                <p>Now published.</p>
              </div>

              <figure className="figure-two-preview">
                <a href="/figures/Figure2.pdf" download="Figure2.pdf">
                  <img src="/figures/Figure2.png" alt="Figure 2 high-resolution preview" />
                </a>
              </figure>

              <p className="figure-caption">
                Fig. 2. Maximum-likelihood (ML) phylogeny of cyanobacterial lineages showing the presence or absence
                of <em>nifHDKENB</em>.
              </p>

              <figure className="figure-two-preview">
                <a href="/figures/Figure3.pdf" download="Figure3.pdf">
                  <img src="/figures/Figure3.png" alt="Figure 3 high-resolution preview" />
                </a>
              </figure>

              <p className="figure-caption">
                Fig. 3. Maximum-likelihood phylogenetic tree of Group I-IV Nif, Vnf, and Anf genes.
              </p>

              <figure className="figure-two-preview">
                <a href="/figures/Figure4.pdf" download="Figure4.pdf">
                  <img src="/figures/Figure4.png" alt="Figure 4 high-resolution preview" />
                </a>
              </figure>

              <p className="figure-caption">
                Fig. 4. Comparison of Groups I and II <em>nif</em> clusters in cyanobacteria.
              </p>

              <figure className="figure-two-preview">
                <a href="/figures/Figure7.pdf" download="Figure7.pdf">
                  <img src="/figures/Figure7.png" alt="Figure 7 high-resolution preview" />
                </a>
              </figure>

              <p className="figure-caption">
                Fig. 7. Metagenomic profiling of cyanobacteria at the genus level across diverse environments.
              </p>

              <figure className="figure-two-preview">
                <a href="/figures/Supplementary_Figure12.pdf" download="Supplementary_Figure12.pdf">
                  <img src="/figures/Supplementary_Figure12.png" alt="Supplementary Figure 12 high-resolution preview" />
                </a>
              </figure>

              <p className="figure-caption">Supplementary Fig. 12.</p>
            </article>
          ) : null}

          {activeTab === "manual" ? (
            <article className="manual-body">
              <h2>Manual</h2>
              <p>
                Nif-Finder detects and classifies nitrogen fixation genes, including <em>nifH</em>, <em>nifD</em>,{" "}
                <em>nifK</em>, <em>nifE</em>, <em>nifN</em>, <em>nifB</em>, <em>vnfH/nifH</em>, <em>vnfD</em>, <em>vnfK</em>, <em>vnfE/nifE</em>, <em>vnfN/nifN</em>, <em>vnfG</em>,{" "}
                <em>vnfDG</em>, <em>nifZ</em>, <em>nifX</em>, <em>nifP/cysE</em>, <em>nifT</em>, <em>nifV</em>, <em>nifS</em>, <em>nifU</em>, <em>nifU_like</em>, <em>modAlike</em>,{" "}
                <em>modB/vupB</em>, <em>modC/vupC</em>, <em>vupA/modA</em>, <em>vupB/modB</em>, <em>vupC/modC</em>, and <em>cnfR/patB</em>, from protein FASTA input. It uses
                HMMER3 hmmscan and nearest-neighbour classification on homology and protein length plots.
              </p>

              <figure className="manual-figure">
                <img src="/manual-run-annotated.jpg" alt="Annotated Nif-Finder Run page controls" />
              </figure>

              <h3>1. Paste or upload protein FASTA</h3>
              <p>
                Paste protein sequences into the Protein FASTA field. Instead of pasting sequences, you can upload a
                local <code>.faa</code>, <code>.fa</code>, <code>.fasta</code>, or <code>.txt</code> file. The web
                interface accepts only protein FASTA input, up to 10 MB. For a demo, you can load an example dataset.
              </p>

              <h3>2. Upload GenBank for genome position plots (optional)</h3>
              <p>
                When a full GenBank file is provided, Nif-Finder shows where the detected target{" "}
                genes are located on the genome. This is useful for checking gene order, fragmented genes, and neighbouring coding
                sequences. The GenBank upload limit is 30 MB.
              </p>

              <h3>3. Adjust analysis parameters (optional)</h3>
              <p>
                Jobs and CPU are parameters for the homology search. The default settings are jobs = 1, CPU = 4, and
                E-value threshold = 1e-10.
              </p>

              <h3>4. Choose output options (optional)</h3>
              <p>
                Plot output controls whether the scatter plot is returned.
              </p>

              <h3>5. Run analysis</h3>
              <p>
                Press Run analysis to start the analysis. A run usually takes 1 to 5 minutes. Results include the query
                identifier, -log10(E-value), alignment length, query protein length, predicted gene, and completeness
                status. Target hits are labelled as Full, Fragment, or Operon. The ZIP download contains TSV and CSV result
                tables, detected target FASTA sequences, the scatter plot, any genome context figures generated from the
                optional GenBank input, and an annotated GenBank file for the visualized local context region when
                available.
              </p>

              <h2 className="manual-section-break">Result</h2>

              <h3>1. The scatter plot of the nif/vnf and related homologues</h3>
              <p>
                The scatter plot summarizes homology search results for target <em>nif</em>/<em>vnf</em> and related proteins. Each panel
                corresponds to one protein: <em>nifH</em>, <em>nifD</em>, <em>nifK</em>, <em>nifE</em>,{" "}
                <em>nifN</em>, <em>nifB</em>, <em>vnfH</em>, <em>vnfD</em>, <em>vnfK</em>, <em>vnfE</em>, <em>vnfN</em>, <em>vnfG</em>, <em>nifZ</em>, <em>nifX</em>, <em>nifP</em>, <em>nifT</em>, <em>nifV</em>, <em>nifS</em>,{" "}
                <em>nifU</em>, <em>modB</em>, <em>modC</em>, <em>modA</em>, <em>vupA</em>, <em>vupB</em>, <em>vupC</em>, or <em>cnfR/patB</em>, with <em>vnfDG</em> shown on the VnfG panel and like labels shown on their parent panels. The x-axis shows protein length in amino acids, and the y-axis shows
                -log10(E-value), so points higher on the plot represent stronger matches. Full-length hits of the target
                proteins are plotted in and around the dashed circle. Circle plots represent hits from
                complete genomes, and triangle plots represent hits from draft genomes. Partial-length points outside
                the full-length region may indicate fragmented <em>nif</em>/<em>vnf</em> genes, especially in draft assemblies.
                Double dashed circles show the <em>nifEN</em> operon.
              </p>
              <figure className="manual-figure manual-figure-wide">
                <img src="/manual-scatter-results.jpg" alt="Annotated Nif-Finder scatter plot result explanation" />
              </figure>

              <h3>2. Visualization of the nif/vnf-related region</h3>
              <p>
                When a GenBank file is provided, Nif-Finder draws genome position plots for detected{" "}
                <em>nifHDKENB</em>, supported <em>vnf</em> genes, and related target genes. The overview plot shows where the target region is encoded on the full
                sequence, and the local plot shows a pinpoint view of the target region. Colored arrows
                show detected target genes, and gray arrows show other ORFs based on the user-provided GenBank
                file. The labels below target genes show whether each hit is full-length, fragmented,
                or part of an operon. In the pinpoint view, positions are recalculated from the left edge of the
                enlarged region, with the upstream end set to +1. These plots are useful for checking probable{" "}
                <em>nif</em>/<em>vnf</em>-cluster regions and nearby target genes for cluster comparison using clinker.
              </p>
              <figure className="manual-figure manual-figure-wide">
                <img src="/manual-genome-context-results.jpg" alt="Nif-Finder genome position and local context result" />
              </figure>
            </article>
          ) : null}

          {activeTab === "about" ? (
            <div className="about-list">
              <article className="info-card citation-card">
                <div>
                  <h3>Citation</h3>
                  <p>
                    Uesaka K, Fujita Y. Accurate prediction of nitrogen fixation in cyanobacteria reveals the dynamic
                    evolution driving high retention rate with mosaic distribution. <em>bioRxiv</em>. 2026.
                  </p>
                  <a href="https://doi.org/10.64898/2026.01.15.699626" target="_blank" rel="noreferrer">
                    https://doi.org/10.64898/2026.01.15.699626
                  </a>
                  <a
                    className="github-icon-link"
                    href="https://github.com/kazumaxneo/Nif_finder"
                    target="_blank"
                    rel="noreferrer"
                    aria-label="Open Nif-Finder on GitHub"
                    title="Open Nif-Finder on GitHub"
                  >
                    <img src="/github-mark.png" alt="" aria-hidden="true" />
                  </a>
                </div>
                </article>
            </div>
          ) : null}
        </section>
      )}
    </main>
  );
}
