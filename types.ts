
export type TargetCYP = 'CYP3A4' | 'CYP2C9' | 'CYP1A2' | 'CYP2D6' | 'CYP2C19' | 'CYP3A5';

export interface GeneInfo {
  symbol: TargetCYP;
  name: string;
  entrez_id: string;
  uniprot_id: string | null;
}

export interface SequenceRecord {
  id: string;
  description: string;
  seq: string;
}

export interface Variant {
  hgvs: string;
  proteinChange: string;
  residue: number;
  targetAA: string;
  refAA: string;
  amScore: number | null; // AlphaMissense Score
  clinVarId?: string;
  clinVarVariantId?: string | number;
  gnomadFreq?: number | null;
  gnomadFreqLabel?: string;
  rsid?: string;
  gnomadLink?: string | null;
  clinicalSignificance?: string;
  distanceToDrug?: number | null;
}

export interface PipelineState {
  step: 'idle' | 'searching' | 'analyzing' | 'designing' | 'complete' | 'error';
  error?: string;
  logs: string[];
}

export interface AdvancedSettings {
  crispr: {
    disruptionPriority: 'PAM' | 'SEED' | 'BOTH';
    seedLength: number;
    minSeedMutations: number;
    maxSeedMutations: number;
    pamConstraint: 'NGG' | 'NNGRRT' | 'TTTV' | 'NG';
    guideLength: number;
    minDoenchScore: number;
    repairTemplateLength: number;
    repairAsymmetry: 'CENTERED' | 'UPSTREAM_SKEW' | 'DOWNSTREAM_SKEW';
    primerSizeMin: number;
    primerSizeMax: number;
    primerTmMin: number;
    primerTmMax: number;
    forceGcClamp: boolean;
  };
  ai: {
    labResources: string[];
    assayPreference: string;
    safetyLevel: 'STANDARD' | 'CLASSROOM_SAFE';
  };
  structure: {
    defaultRepresentation: 'cartoon' | 'stick' | 'surface';
    colorScheme: 'chain' | 'secondary' | 'conservation';
  };
}

// --- CRISPR Types ---

export interface Cas9Site {
  position: number;
  sequence: string; // The full match including PAM
  strand: 'forward' | 'reverse';
  context30?: string; // 4bp + 20bp + PAM + 3bp for scoring
}

export interface VerificationPrimers {
  forward: string;
  reverse: string;
  productSize: number;
  forwardTm: number;
  reverseTm: number;
  forwardStart: number;
  reverseStart: number;
}

export interface RepairResult {
  site: Cas9Site;
  cloningOligoA: string;
  cloningOligoB: string;
  repairTemplate: string; // The mutated homology arm
  guideSeqWithPam: string; // N20NGG format (5'->3')
  score?: number; // Efficiency Score (0-100)
  
  // Deletion Control Fields
  deletionRepairTemplate: string;
  deletionDnaDisplay: string; // Original sequence with dashes where PAM was
  deletionProtein: string; // Translated protein of the deletion mutant

  originalRegion: string;
  homologyStart: number;
  mutationPosition: number; // Nucleotide index relative to start of gene
  aaChangeStatus: 'success' | 'warning';
  aaChangesCount: number;
  dnaAlignment: AlignmentData;
  aaAlignment: AlignmentData;
  strategy: 'PAM_DISRUPTED_BY_TARGET' | 'PAM_SILENT' | 'SEED_SILENT';
  silentMutationCount: number;
  
  verificationPrimers?: VerificationPrimers;
}

export interface AlignmentData {
  original: string;
  modified: string;
  matchString: string; // String of '|' and ' '
}

export interface CodonTable {
  [key: string]: string[];
}

export interface CodingExon {
  start: number; // String index (inclusive)
  end: number;   // String index (exclusive)
  cumLengthBefore: number; // Cumulative coding bases before this exon
}

