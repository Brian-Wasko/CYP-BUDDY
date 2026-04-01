
import { GeneInfo, SequenceRecord, TargetCYP } from '../types';

// --- MyGene.info ---
export const getHumanGeneInfo = async (symbol: TargetCYP): Promise<GeneInfo> => {
  const url = `https://mygene.info/v3/query?q=${symbol}&scopes=symbol,alias&fields=symbol,entrezgene,name,uniprot,type_of_gene&species=human`;
  
  const response = await fetch(url);
  const data = await response.json();
  
  if (!data.hits || data.hits.length === 0) {
    throw new Error("Gene not found in MyGene.info");
  }
  const hit = data.hits[0];

  let uniprotId = null;
  if (hit.uniprot) {
    if (typeof hit.uniprot === 'string') uniprotId = hit.uniprot;
    else if (hit.uniprot['Swiss-Prot']) uniprotId = Array.isArray(hit.uniprot['Swiss-Prot']) ? hit.uniprot['Swiss-Prot'][0] : hit.uniprot['Swiss-Prot'];
  }

  return {
    symbol: hit.symbol as TargetCYP,
    name: hit.name,
    entrez_id: hit.entrezgene?.toString(),
    uniprot_id: uniprotId
  };
};

// --- UniProt ---
export const fetchSequence = async (id: string): Promise<SequenceRecord> => {
  const url = `https://rest.uniprot.org/uniprotkb/${id}.fasta`;
  
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`Failed to fetch sequence for ${id} (Status: ${response.status})`);
  }
  
  const text = await response.text();
  if (!text || !text.trim()) {
     throw new Error(`No sequence data returned for ${id}`);
  }

  const lines = text.trim().split('\n');
  const description = lines[0];
  const seq = lines.slice(1).join('').trim();
  
  return { id, description, seq };
};

// --- Ensembl ---
const ENSEMBL_IDS: Record<TargetCYP, string> = {
  'CYP3A4': 'ENSG00000160868',
  'CYP2C9': 'ENSG00000138109',
  'CYP1A2': 'ENSG00000140505',
  'CYP2D6': 'ENSG00000100197',
  'CYP2C19': 'ENSG00000165841',
  'CYP3A5': 'ENSG00000106258'
};

export const fetchHumanCdna = async (symbol: TargetCYP): Promise<string> => {
  const ensemblId = ENSEMBL_IDS[symbol];
  if (!ensemblId) throw new Error(`Ensembl ID not found for ${symbol}`);

  // Fetch transcripts
  const url = `https://rest.ensembl.org/lookup/id/${ensemblId}?expand=1`;
  const response = await fetch(url, { headers: { 'Content-Type': 'application/json' } });
  if (!response.ok) throw new Error(`Failed to fetch Ensembl data for ${symbol}`);
  
  const data = await response.json();
  const transcripts = data.Transcript || [];
  
  // Find canonical or first protein-coding
  let targetTranscript = transcripts.find((t: any) => t.is_canonical) || transcripts.find((t: any) => t.biotype === 'protein_coding') || transcripts[0];
  
  if (!targetTranscript) throw new Error(`No transcript found for ${symbol}`);

  // Fetch CDS sequence
  const seqUrl = `https://rest.ensembl.org/sequence/id/${targetTranscript.id}?type=cds`;
  const seqResponse = await fetch(seqUrl, { headers: { 'Content-Type': 'text/plain' } });
  if (!seqResponse.ok) throw new Error(`Failed to fetch CDS for ${targetTranscript.id}`);
  
  const cdna = await seqResponse.text();
  return cdna.trim();
};
export const fetchClinVarVariants = async (geneSymbol: string, significanceTerms: string[] = []): Promise<any[]> => {
  const sigQuery = significanceTerms.length > 0 
      ? `(${significanceTerms.map(s => `clinvar.rcv.clinical_significance:"${s}"`).join(" OR ")})`
      : `clinvar.rcv.clinical_significance:"Uncertain significance"`;

  const query = `clinvar.gene.symbol:${geneSymbol} AND ${sigQuery}`;
  const url = `https://myvariant.info/v1/query?q=${query}&fields=clinvar,dbsnp,dbnsfp,gnomad_exome,gnomad_genome&size=1000`;
  
  const response = await fetch(url);
  const data = await response.json();
  return data.hits || [];
};

