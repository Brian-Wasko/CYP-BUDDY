import React, { useState, useEffect, useRef, useMemo } from 'react';
import { Search, Dna, Activity, FileText, AlertCircle, PlayCircle, Settings, ExternalLink, Info, List, ArrowRight, Sparkles, Filter, FlaskConical, Copy, Download, HelpCircle, ChevronDown, Mail, X, Printer, Loader2, ChevronLeft, ChevronRight, Sun, Moon, RefreshCw, Link2, Box, ArrowUp, ArrowDown } from 'lucide-react';
import ReactMarkdown from 'react-markdown';
import { GeneInfo, Variant, PipelineState, RepairResult, Cas9Site, AdvancedSettings, TargetCYP, SequenceRecord } from './types';
import { getHumanGeneInfo, fetchSequence, fetchClinVarVariants, fetchHumanCdna } from './services/api';
import { parseProteinChange } from './utils/alignment';
import { findCas9Sites, generateRepairTemplates, getMutationIndex } from './utils/crispr';
import { fetchPdbForGeneAndDrug, calculateDistancesToDrug } from './utils/pdb';
import { generateExperimentalPlan, fetchSuggestedDrugs } from './services/geminiService';
import { StructureViewer, StructureViewerHandle } from './components/StructureViewer';

const CYP_GENES: TargetCYP[] = ['CYP3A4', 'CYP2C9', 'CYP1A2', 'CYP2D6', 'CYP2C19', 'CYP3A5'];

const DEFAULT_SETTINGS: AdvancedSettings = {
    crispr: {
        disruptionPriority: 'BOTH',
        seedLength: 10,
        minSeedMutations: 2,
        maxSeedMutations: 5,
        pamConstraint: 'NGG',
        guideLength: 20,
        minDoenchScore: 0,
        repairTemplateLength: 80,
        repairAsymmetry: 'CENTERED',
        primerSizeMin: 450,
        primerSizeMax: 900,
        primerTmMin: 53,
        primerTmMax: 62,
        forceGcClamp: false
    },
    ai: {
        labResources: [],
        assayPreference: 'ANY',
        safetyLevel: 'CLASSROOM_SAFE'
    },
    structure: {
        defaultRepresentation: 'cartoon',
        colorScheme: 'chain'
    }
};

const renderRefDna = (seq: string, site: Cas9Site, homologyStart: number) => {
  const pamStartRel = site.position - homologyStart;
  let pamStart = -1;
  let pamEnd = -1;

  if (site.strand === 'forward') {
    pamStart = pamStartRel + 20;
    pamEnd = pamStart + 3;
  } else {
    pamStart = pamStartRel;
    pamEnd = pamStart + 3;
  }

  return (
    <span>
      {seq.split('').map((char, i) => {
         const isPam = i >= pamStart && i < pamEnd;
         return isPam ? <span key={i} className="text-purple-600 dark:text-purple-400 font-bold">{char}</span> : char;
      })}
    </span>
  );
};

const renderAltSeq = (ref: string, alt: string) => {
   return (
     <span>
       {alt.split('').map((char, i) => {
         const refChar = ref[i] || '';
         const isDiff = char !== refChar;
         return isDiff ? <span key={i} className="text-red-600 dark:text-red-400 font-bold">{char}</span> : char;
       })}
     </span>
   );
};

const getFreqLabel = (af: number | null | undefined) => {
  if (af === null || af === undefined) return 'Unknown';
  if (af >= 0.05) return 'Common';
  if (af >= 0.01) return 'Low Freq';
  if (af >= 0.001) return 'Rare';
  return 'Very Rare';
};

const parseAmScore = (score: any): number | null => {
  if (score === null || score === undefined) return null;
  if (Array.isArray(score)) {
    return score.length > 0 ? Number(score[0]) : null;
  }
  if (typeof score === 'string') {
    const parts = score.split(',');
    return Number(parts[0]);
  }
  const num = Number(score);
  return isNaN(num) ? null : num;
};

export const App: React.FC = () => {
  const [isDarkMode, setIsDarkMode] = useState(true);
  const [selectedGene, setSelectedGene] = useState<TargetCYP>('CYP3A4');
  const [drugName, setDrugName] = useState('');
  const [variantInput, setVariantInput] = useState('');
  
  const [suggestedDrugs, setSuggestedDrugs] = useState<string[]>([]);
  const [suggestedVariants, setSuggestedVariants] = useState<Variant[]>([]);
  const [isSearchingSuggestions, setIsSearchingSuggestions] = useState(false);

  type SortColumn = 'proteinChange' | 'clinicalSignificance' | 'amScore' | 'distanceToDrug' | 'gnomadFreq';
  const [sortColumn, setSortColumn] = useState<SortColumn>('amScore');
  const [sortDirection, setSortDirection] = useState<'asc' | 'desc'>('desc');

  const sortedVariants = useMemo(() => {
    return [...suggestedVariants].sort((a, b) => {
      let valA: any = a[sortColumn];
      let valB: any = b[sortColumn];

      if (sortColumn === 'amScore' || sortColumn === 'distanceToDrug' || sortColumn === 'gnomadFreq') {
        valA = valA !== null && valA !== undefined && !isNaN(Number(valA)) ? Number(valA) : (sortDirection === 'asc' ? Infinity : -Infinity);
        valB = valB !== null && valB !== undefined && !isNaN(Number(valB)) ? Number(valB) : (sortDirection === 'asc' ? Infinity : -Infinity);
      } else {
        valA = String(valA || '').toLowerCase();
        valB = String(valB || '').toLowerCase();
      }

      if (valA < valB) return sortDirection === 'asc' ? -1 : 1;
      if (valA > valB) return sortDirection === 'asc' ? 1 : -1;
      return 0;
    });
  }, [suggestedVariants, sortColumn, sortDirection]);

  const handleSort = (column: SortColumn) => {
    if (sortColumn === column) {
      setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
    } else {
      setSortColumn(column);
      setSortDirection('desc'); // Default to desc for scores, asc for text
    }
  };

  const [state, setState] = useState<PipelineState>({ step: 'idle', logs: [] });
  const [isGeneratingPlan, setIsGeneratingPlan] = useState(false);
  const [isGeneratingCrispr, setIsGeneratingCrispr] = useState(false);
  const [geneInfo, setGeneInfo] = useState<GeneInfo | null>(null);
  const [sequence, setSequence] = useState<any>(null);
  const [cdna, setCdna] = useState<string>('');
  const [variant, setVariant] = useState<Variant | null>(null);
  const [aiPlan, setAiPlan] = useState<string>('');
  const [crisprResults, setCrisprResults] = useState<RepairResult[]>([]);
  const [selectedCrisprIndex, setSelectedCrisprIndex] = useState(0);
  const [isGeneratingDiffDock, setIsGeneratingDiffDock] = useState(false);
  const [diffDockResults, setDiffDockResults] = useState<any>(null);

  useEffect(() => {
    setVariant(prev => prev && prev.distanceToDrug !== undefined ? { ...prev, distanceToDrug: undefined } : prev);
    setSuggestedVariants(prev => {
      if (prev.some(v => v.distanceToDrug !== undefined)) {
        return prev.map(v => ({ ...v, distanceToDrug: undefined }));
      }
      return prev;
    });
  }, [drugName]);

  const structureViewerRef = useRef<StructureViewerHandle>(null);
  const logsEndRef = useRef<HTMLDivElement>(null);

  const fetchSuggestions = async () => {
    setIsSearchingSuggestions(true);
    setSuggestedDrugs([]);
    setSuggestedVariants([]);
    
    try {
      // 1. Fetch Drugs via Gemini
      const drugs = await fetchSuggestedDrugs(selectedGene);
      setSuggestedDrugs(drugs);

      // 2. Fetch Variants via MyVariant.info
      const rawHits = await fetchClinVarVariants(selectedGene, ['Uncertain significance', 'Pathogenic', 'Likely pathogenic']);
      let variants: Variant[] = rawHits.slice(0, 15).map(hit => {
        let clinVarEntry = hit.clinvar;
        if (Array.isArray(clinVarEntry)) clinVarEntry = clinVarEntry[0];
        const pChange = clinVarEntry?.hgvs?.protein;
        const pChangeStr = Array.isArray(pChange) ? pChange[0] : pChange;
        
        const parsed = pChangeStr ? parseProteinChange(pChangeStr) : null;
        const shortChangeStr = pChangeStr ? (pChangeStr.split('p.')[1] || pChangeStr) : 'N/A';
        
        let gnomadFreq = null;
        if (hit.gnomad_exome?.af?.af) gnomadFreq = hit.gnomad_exome.af.af;
        else if (hit.gnomad_genome?.af?.af) gnomadFreq = hit.gnomad_genome.af.af;

        const rsid = hit.dbsnp?.rsid;
        const variantId = hit._id; // usually chr:pos:ref:alt or similar

        return {
          hgvs: pChangeStr || 'N/A',
          proteinChange: shortChangeStr,
          residue: parsed?.res || 0,
          targetAA: parsed?.target || '?',
          refAA: parsed?.ref || '?',
          amScore: parseAmScore(hit.dbnsfp?.alphamissense?.score),
          clinicalSignificance: clinVarEntry?.rcv?.clinical_significance || 'Unknown',
          clinVarId: clinVarEntry?.rcv?.accession,
          gnomadFreq: gnomadFreq,
          gnomadFreqLabel: getFreqLabel(gnomadFreq),
          rsid: rsid,
          gnomadLink: variantId ? `https://gnomad.broadinstitute.org/search?query=${rsid || variantId}` : null
        };
      }).filter(v => v.proteinChange !== 'N/A');
      
      // 3. Calculate distances if drug is present
      const targetDrug = drugName.trim() || drugs[0];
      if (targetDrug) {
        try {
          const gInfo = await getHumanGeneInfo(selectedGene);
          if (gInfo.uniprot_id) {
            const pdbData = await fetchPdbForGeneAndDrug(gInfo.uniprot_id, targetDrug);
            if (pdbData) {
              const residues = variants.map(v => v.residue).filter(r => r > 0);
              const distances = calculateDistancesToDrug(pdbData, residues);
              variants = variants.map(v => ({
                ...v,
                distanceToDrug: distances[v.residue] !== undefined ? distances[v.residue] : null
              }));
            }
          }
        } catch (e) {
          console.warn("Could not calculate distances", e);
        }
      }

      setSuggestedVariants(variants);
    } catch (error) {
      console.error("Error fetching suggestions:", error);
    } finally {
      setIsSearchingSuggestions(false);
    }
  };

  useEffect(() => {
    if (isDarkMode) document.documentElement.classList.add('dark');
    else document.documentElement.classList.remove('dark');
  }, [isDarkMode]);

  useEffect(() => {
    logsEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [state.logs, state.error]);

  const addLog = (msg: string) => {
    setState(prev => ({ ...prev, logs: [...prev.logs, msg] }));
  };

  const handleAnalyze = async () => {
    if (!drugName.trim() || !variantInput.trim()) {
      if (suggestedDrugs.length === 0 && suggestedVariants.length === 0) {
        await fetchSuggestions();
        setState({ step: 'idle', logs: ["Please select a drug and variant from the suggestions below or enter them manually."] });
      } else {
        setState({ step: 'error', error: 'Please enter both a drug name and a variant.', logs: [] });
      }
      return;
    }

    setState({ step: 'searching', logs: ['Starting analysis pipeline...'] });
    setGeneInfo(null);
    setVariant(null);
    setAiPlan('');
    setCrisprResults([]);

    try {
      // 1. Parse Variant
      const parsed = parseProteinChange(variantInput);
      if (!parsed) {
        throw new Error("Invalid variant format. Please use format like 'R144C' or 'p.Arg144Cys'.");
      }

      // 2. Fetch Gene Info
      addLog(`Fetching human gene info for ${selectedGene}...`);
      const gInfo = await getHumanGeneInfo(selectedGene);
      setGeneInfo(gInfo);

      // 3. Fetch Sequence & cDNA
      addLog(`Fetching sequences for ${selectedGene}...`);
      const seqRecord = await fetchSequence(gInfo.uniprot_id!);
      setSequence(seqRecord);
      const cdsSeq = await fetchHumanCdna(selectedGene);
      setCdna(cdsSeq);

      // 4. Fetch ClinVar Data
      addLog(`Checking ClinVar for ${selectedGene} variants...`);
      const rawHits = await fetchClinVarVariants(selectedGene, []);
      let foundVariant: Variant = {
        hgvs: variantInput,
        proteinChange: variantInput.split('p.')[1] || variantInput,
        residue: parsed.res,
        targetAA: parsed.target,
        refAA: parsed.ref,
        amScore: null,
        clinicalSignificance: 'Unknown',
        gnomadFreq: null
      };

      for (const hit of rawHits) {
        let clinVarEntry = hit.clinvar;
        if (Array.isArray(clinVarEntry)) clinVarEntry = clinVarEntry[0];
        const pChange = clinVarEntry?.hgvs?.protein;
        const pChangeStr = Array.isArray(pChange) ? pChange[0] : pChange;

        if (pChangeStr) {
          const hitParsed = parseProteinChange(pChangeStr);
          if (hitParsed && hitParsed.res === parsed.res && hitParsed.target === parsed.target) {
             let gnomadFreq = null;
             if (hit.gnomad_exome?.af?.af) {
                 gnomadFreq = hit.gnomad_exome.af.af;
             } else if (hit.gnomad_genome?.af?.af) {
                 gnomadFreq = hit.gnomad_genome.af.af;
             }

             const rsid = hit.dbsnp?.rsid;
             const variantId = hit._id;

             foundVariant = {
               ...foundVariant,
               hgvs: pChangeStr,
               clinVarId: clinVarEntry.rcv?.accession,
               clinicalSignificance: clinVarEntry.rcv?.clinical_significance,
               amScore: parseAmScore(hit.dbnsfp?.alphamissense?.score),
               gnomadFreq: gnomadFreq,
               gnomadFreqLabel: getFreqLabel(gnomadFreq),
               rsid: rsid,
               gnomadLink: variantId ? `https://gnomad.broadinstitute.org/search?query=${rsid || variantId}` : null
             };
             break;
          }
        }
      }
      
      // 5. Calculate distance if drug is present
      if (drugName.trim() && gInfo.uniprot_id) {
        try {
          addLog(`Checking structure distance for ${drugName}...`);
          const pdbData = await fetchPdbForGeneAndDrug(gInfo.uniprot_id, drugName.trim());
          if (pdbData) {
            const distances = calculateDistancesToDrug(pdbData, [foundVariant.residue]);
            if (distances[foundVariant.residue] !== undefined) {
              foundVariant.distanceToDrug = distances[foundVariant.residue];
              addLog(`Distance to drug calculated: ${foundVariant.distanceToDrug.toFixed(1)} Å`);
            }
          }
        } catch (e) {
          console.warn("Could not calculate distance in handleAnalyze", e);
        }
      }

      setVariant(foundVariant);
      addLog(`Variant processed: ${foundVariant.proteinChange}`);

      setState({ step: 'complete', logs: ['Target loaded. Ready for analysis.'] });

    } catch (err: any) {
      console.error(err);
      setState({ step: 'error', error: err.message, logs: [...state.logs, `Error: ${err.message}`] });
    }
  };

  const handleCalculateDistances = async () => {
    if (!drugName.trim() || !geneInfo?.uniprot_id) {
      addLog("Please select a drug and load a target first.");
      return;
    }
    
    addLog(`Calculating structural distances for ${drugName}...`);
    try {
      const pdbData = await fetchPdbForGeneAndDrug(geneInfo.uniprot_id, drugName.trim());
      if (pdbData) {
        // Calculate for main variant
        if (variant) {
          const distances = calculateDistancesToDrug(pdbData, [variant.residue]);
          if (distances[variant.residue] !== undefined) {
            setVariant({ ...variant, distanceToDrug: distances[variant.residue] });
            addLog(`Distance to drug calculated for ${variant.proteinChange}: ${distances[variant.residue].toFixed(1)} Å`);
          } else {
            addLog(`Could not find distance for ${variant.proteinChange}.`);
          }
        }
        
        // Calculate for suggested variants
        if (suggestedVariants.length > 0) {
          const residues = suggestedVariants.map(v => v.residue);
          const distances = calculateDistancesToDrug(pdbData, residues);
          const updatedSuggestions = suggestedVariants.map(v => ({
            ...v,
            distanceToDrug: distances[v.residue] !== undefined ? distances[v.residue] : v.distanceToDrug
          }));
          setSuggestedVariants(updatedSuggestions);
          addLog(`Distances calculated for suggested variants.`);
        }
      } else {
        addLog(`Could not find co-complex structure for ${selectedGene} and ${drugName}.`);
      }
    } catch (e) {
      console.warn("Could not calculate distance", e);
      addLog(`Error calculating distance: ${e}`);
    }
  };

  const handleGenerateCrispr = () => {
    if (!variant || !cdna) return;
    setIsGeneratingCrispr(true);
    addLog(`Designing CRISPR oligos...`);
    
    try {
      const codingExons = [{ start: 0, end: cdna.length, cumLengthBefore: 0 }];
      const mutIndex = getMutationIndex(codingExons, variant.residue);
      
      const cSites = findCas9Sites(cdna, mutIndex, DEFAULT_SETTINGS.crispr);
      if (cSites.length === 0) {
          addLog("Warning: No suitable Cas9 sites found near the mutation.");
      } else {
          const tmpls = generateRepairTemplates(cdna, cSites, mutIndex, variant.targetAA, DEFAULT_SETTINGS.crispr.repairTemplateLength, DEFAULT_SETTINGS.crispr);
          setCrisprResults(tmpls);
          addLog(`Generated ${tmpls.length} CRISPR repair strategies.`);
      }
    } catch (e) {
      console.error(e);
      addLog("Error generating CRISPR designs.");
    } finally {
      setIsGeneratingCrispr(false);
    }
  };

  const handleGenerateDiffDock = async () => {
    setIsGeneratingDiffDock(true);
    addLog("Making request to DiffDock via backend...");
    try {
      const response = await fetch("/api/diffdock", {
        method: "POST",
        headers: {
          "Content-Type": "application/json"
        }
      });

      if (!response.ok) {
        throw new Error(`DiffDock API error: ${response.statusText}`);
      }

      const data = await response.json();
      setDiffDockResults(data);
      addLog("DiffDock modelling complete.");
    } catch (error: any) {
      console.error(error);
      addLog(`ERROR: ${error.message}`);
    } finally {
      setIsGeneratingDiffDock(false);
    }
  };

  const handleGeneratePlan = async () => {
    if (!variant || !geneInfo) return;
    setIsGeneratingPlan(true);
    addLog(`Generating experimental plan with Gemini...`);
    
    try {
      let structureImg = null;
      if (structureViewerRef.current) {
          structureImg = structureViewerRef.current.captureImage();
      }

      await generateExperimentalPlan(
        selectedGene,
        drugName,
        variant,
        (text) => setAiPlan(text),
        structureImg
      );
      addLog(`Experimental plan generated successfully.`);
    } catch (e) {
      console.error(e);
      addLog("Error generating experimental plan.");
    } finally {
      setIsGeneratingPlan(false);
    }
  };

  const handleVariantSelect = async (v: Variant) => {
    setVariantInput(v.proteinChange);
    setVariant(v);
    setAiPlan('');
    setCrisprResults([]);
    
    // Ensure gene info is loaded if they click a variant before "Set Target"
    if (!geneInfo || geneInfo.symbol !== selectedGene) {
      setState({ step: 'searching', logs: [...state.logs, `Loading target data for ${selectedGene}...`] });
      try {
        const gInfo = await getHumanGeneInfo(selectedGene);
        setGeneInfo(gInfo);
        const seqRecord = await fetchSequence(gInfo.uniprot_id!);
        setSequence(seqRecord);
        const cdsSeq = await fetchHumanCdna(selectedGene);
        setCdna(cdsSeq);
        setState({ step: 'complete', logs: [...state.logs, `Target data loaded for ${selectedGene}. Variant selected: ${v.proteinChange}`] });
      } catch (e) {
        console.error("Failed to load gene info for selected variant", e);
        setState({ step: 'error', error: 'Failed to load target data.', logs: [...state.logs, `Error: Failed to load target data.`] });
      }
    } else {
      addLog(`Variant selected: ${v.proteinChange}`);
      setState(prev => ({ ...prev, step: 'complete' }));
    }
  };

  return (
    <div className="min-h-screen bg-slate-50 dark:bg-slate-900 text-slate-900 dark:text-slate-100 font-sans flex flex-col">
      {/* Header */}
      <header className="bg-white dark:bg-slate-800 border-b border-slate-200 dark:border-slate-700 p-4 flex justify-between items-center z-10 sticky top-0 shadow-sm">
        <div className="flex items-center space-x-3">
          <div className="bg-blue-600 p-2 rounded-lg shadow-inner">
            <FlaskConical className="w-6 h-6 text-white" />
          </div>
          <div>
            <h1 className="text-xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-blue-600 to-purple-600 dark:from-blue-400 dark:to-purple-400">
              CYP Variant Analysis & 12Δ Yeast Assay Designer
            </h1>
            <p className="text-xs text-slate-500 dark:text-slate-400 font-medium tracking-wide">
              NON-ANIMAL MODEL TOXICOLOGY PLATFORM
            </p>
          </div>
        </div>
        <button onClick={() => setIsDarkMode(!isDarkMode)} className="p-2 rounded-full hover:bg-slate-100 dark:hover:bg-slate-700 transition-colors">
          {isDarkMode ? <Sun className="w-5 h-5 text-amber-400" /> : <Moon className="w-5 h-5 text-slate-500" />}
        </button>
      </header>

      <main className="flex-1 max-w-7xl w-full mx-auto p-4 md:p-6 grid grid-cols-1 lg:grid-cols-12 gap-6">
        
        {/* Left Panel: Inputs */}
        <div className="lg:col-span-4 space-y-6">
          <div className="bg-white dark:bg-slate-800 rounded-xl shadow-sm border border-slate-200 dark:border-slate-700 p-5">
            <h2 className="text-lg font-semibold mb-4 flex items-center"><Search className="w-5 h-5 mr-2 text-blue-500" /> Analysis Parameters</h2>
            
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-slate-700 dark:text-slate-300 mb-1">Target Human CYP Gene</label>
                <select 
                  value={selectedGene} 
                  onChange={(e) => {
                    setSelectedGene(e.target.value as TargetCYP);
                    setSuggestedDrugs([]);
                    setSuggestedVariants([]);
                  }}
                  className="w-full bg-slate-50 dark:bg-slate-900 border border-slate-300 dark:border-slate-600 rounded-lg px-4 py-2.5 focus:ring-2 focus:ring-blue-500 focus:border-blue-500 outline-none transition-all"
                >
                  {CYP_GENES.map(g => <option key={g} value={g}>{g}</option>)}
                </select>
              </div>

              <div>
                <label className="block text-sm font-medium text-slate-700 dark:text-slate-300 mb-1">Target FDA-Approved Drug</label>
                <input 
                  type="text" 
                  value={drugName} 
                  onChange={(e) => setDrugName(e.target.value)}
                  placeholder="e.g., Warfarin, Omeprazole"
                  className="w-full bg-slate-50 dark:bg-slate-900 border border-slate-300 dark:border-slate-600 rounded-lg px-4 py-2.5 focus:ring-2 focus:ring-blue-500 focus:border-blue-500 outline-none transition-all"
                />
              </div>

              <div>
                <label className="block text-sm font-medium text-slate-700 dark:text-slate-300 mb-1">Variant (VUS/SNP)</label>
                <input 
                  type="text" 
                  value={variantInput} 
                  onChange={(e) => setVariantInput(e.target.value)}
                  placeholder="e.g., R144C or p.Arg144Cys"
                  className="w-full bg-slate-50 dark:bg-slate-900 border border-slate-300 dark:border-slate-600 rounded-lg px-4 py-2.5 focus:ring-2 focus:ring-blue-500 focus:border-blue-500 outline-none transition-all"
                />
              </div>

              <div className="flex flex-col space-y-2">
                <button 
                  onClick={fetchSuggestions}
                  disabled={isSearchingSuggestions}
                  className="w-full bg-slate-100 dark:bg-slate-700 hover:bg-slate-200 dark:hover:bg-slate-600 text-slate-700 dark:text-slate-200 font-medium py-2 px-4 rounded-lg transition-colors flex items-center justify-center disabled:opacity-50 shadow-sm text-sm"
                >
                  {isSearchingSuggestions ? (
                    <><Loader2 className="w-4 h-4 mr-2 animate-spin" /> Searching Suggestions...</>
                  ) : (
                    <><Sparkles className="w-4 h-4 mr-2 text-amber-500" /> Suggest Drugs & Variants</>
                  )}
                </button>

                {suggestedDrugs.length > 0 && (
                  <div className="bg-slate-50 dark:bg-slate-900/50 p-3 rounded-lg border border-slate-200 dark:border-slate-700">
                    <span className="text-[10px] font-bold text-slate-500 uppercase tracking-wider block mb-2">Suggested Drugs</span>
                    <div className="flex flex-wrap gap-1.5">
                      {suggestedDrugs.map(d => (
                        <button 
                          key={d} 
                          onClick={() => setDrugName(d)}
                          className={`text-[10px] px-2 py-1 rounded transition-all ${
                            drugName === d 
                              ? 'bg-blue-100 dark:bg-blue-900/50 border border-blue-400 dark:border-blue-500 text-blue-700 dark:text-blue-300 font-bold' 
                              : 'bg-white dark:bg-slate-800 border border-slate-300 dark:border-slate-600 hover:bg-blue-50 dark:hover:bg-blue-900/30 hover:border-blue-300'
                          }`}
                        >
                          {d}
                        </button>
                      ))}
                    </div>
                  </div>
                )}

                {suggestedVariants.length > 0 && (
                  <div className="bg-slate-50 dark:bg-slate-900/50 p-3 rounded-lg border border-slate-200 dark:border-slate-700">
                    <div className="flex justify-between items-center mb-2">
                      <span className="text-[10px] font-bold text-slate-500 uppercase tracking-wider block">Suggested Variants</span>
                      {drugName && geneInfo?.uniprot_id && (
                        <button 
                          onClick={handleCalculateDistances}
                          className="text-[10px] bg-blue-100 dark:bg-blue-900/40 text-blue-700 dark:text-blue-300 px-2 py-0.5 rounded hover:bg-blue-200 dark:hover:bg-blue-800/60 transition-colors"
                          title="Calculate distance from variants to drug"
                        >
                          Calc Distances
                        </button>
                      )}
                    </div>
                    <div className="overflow-x-auto">
                      <table className="w-full text-left text-[10px] text-slate-600 dark:text-slate-300">
                        <thead className="bg-slate-100 dark:bg-slate-800 border-b border-slate-200 dark:border-slate-700">
                          <tr>
                            <th className="p-1.5 font-semibold cursor-pointer hover:bg-slate-200 dark:hover:bg-slate-700 transition-colors" onClick={() => handleSort('proteinChange')}>
                              <div className="flex items-center gap-1">Variant {sortColumn === 'proteinChange' && (sortDirection === 'asc' ? <ArrowUp className="w-3 h-3" /> : <ArrowDown className="w-3 h-3" />)}</div>
                            </th>
                            <th className="p-1.5 font-semibold cursor-pointer hover:bg-slate-200 dark:hover:bg-slate-700 transition-colors" onClick={() => handleSort('clinicalSignificance')}>
                              <div className="flex items-center gap-1">Type {sortColumn === 'clinicalSignificance' && (sortDirection === 'asc' ? <ArrowUp className="w-3 h-3" /> : <ArrowDown className="w-3 h-3" />)}</div>
                            </th>
                            <th className="p-1.5 font-semibold cursor-pointer hover:bg-slate-200 dark:hover:bg-slate-700 transition-colors" onClick={() => handleSort('amScore')}>
                              <div className="flex items-center gap-1">AM Score {sortColumn === 'amScore' && (sortDirection === 'asc' ? <ArrowUp className="w-3 h-3" /> : <ArrowDown className="w-3 h-3" />)}</div>
                            </th>
                            <th className="p-1.5 font-semibold cursor-pointer hover:bg-slate-200 dark:hover:bg-slate-700 transition-colors" title="Distance to Drug" onClick={() => handleSort('distanceToDrug')}>
                              <div className="flex items-center gap-1">Dist (Å) {sortColumn === 'distanceToDrug' && (sortDirection === 'asc' ? <ArrowUp className="w-3 h-3" /> : <ArrowDown className="w-3 h-3" />)}</div>
                            </th>
                            <th className="p-1.5 font-semibold cursor-pointer hover:bg-slate-200 dark:hover:bg-slate-700 transition-colors" onClick={() => handleSort('gnomadFreq')}>
                              <div className="flex items-center gap-1">Freq {sortColumn === 'gnomadFreq' && (sortDirection === 'asc' ? <ArrowUp className="w-3 h-3" /> : <ArrowDown className="w-3 h-3" />)}</div>
                            </th>
                            <th className="p-1.5 font-semibold">Links</th>
                          </tr>
                        </thead>
                        <tbody>
                          {sortedVariants.map(v => (
                            <tr 
                              key={v.proteinChange} 
                              onClick={() => handleVariantSelect(v)}
                              className={`border-b border-slate-100 dark:border-slate-800 hover:bg-blue-50 dark:hover:bg-blue-900/20 cursor-pointer transition-colors ${variant?.proteinChange === v.proteinChange ? 'bg-blue-50 dark:bg-blue-900/30' : ''}`}
                            >
                              <td className="p-1.5 font-bold text-blue-600 dark:text-blue-400">{v.proteinChange}</td>
                              <td className="p-1.5 truncate max-w-[80px]" title={v.clinicalSignificance}>{v.clinicalSignificance}</td>
                              <td className="p-1.5">{v.amScore !== null && v.amScore !== undefined && !isNaN(Number(v.amScore)) ? Number(v.amScore).toFixed(2) : '-'}</td>
                              <td className="p-1.5 font-mono">{v.distanceToDrug !== null && v.distanceToDrug !== undefined ? v.distanceToDrug.toFixed(1) : '-'}</td>
                              <td className="p-1.5">
                                {v.gnomadFreq !== null && v.gnomadFreq !== undefined && !isNaN(Number(v.gnomadFreq)) ? (
                                  <span title={Number(v.gnomadFreq).toExponential(2)}>{v.gnomadFreqLabel}</span>
                                ) : '-'}
                              </td>
                              <td className="p-1.5 flex gap-1">
                                {v.clinVarId && <a href={`https://www.ncbi.nlm.nih.gov/clinvar/variation/${v.clinVarId}/`} target="_blank" rel="noreferrer" className="text-blue-500 hover:underline" onClick={e => e.stopPropagation()}>CV</a>}
                                {v.gnomadLink && <a href={v.gnomadLink} target="_blank" rel="noreferrer" className="text-blue-500 hover:underline" onClick={e => e.stopPropagation()}>GN</a>}
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </div>
                )}

                {!isSearchingSuggestions && suggestedDrugs.length === 0 && suggestedVariants.length === 0 && state.logs.includes("Please select a drug and variant from the suggestions below or enter them manually.") && (
                  <div className="text-xs text-amber-600 dark:text-amber-400 bg-amber-50 dark:bg-amber-900/20 p-2 rounded border border-amber-200 dark:border-amber-800 flex items-center">
                    <AlertCircle className="w-3 h-3 mr-2" />
                    No common drugs or variants found for this gene.
                  </div>
                )}
              </div>

              <button 
                onClick={handleAnalyze}
                disabled={state.step !== 'idle' && state.step !== 'complete' && state.step !== 'error'}
                className="w-full bg-blue-600 hover:bg-blue-700 text-white font-medium py-3 px-4 rounded-lg transition-colors flex items-center justify-center disabled:opacity-50 disabled:cursor-not-allowed shadow-sm"
              >
                {state.step === 'searching' ? (
                  <><Loader2 className="w-5 h-5 mr-2 animate-spin" /> Loading Target...</>
                ) : (
                  <><PlayCircle className="w-5 h-5 mr-2" /> Load Target</>
                )}
              </button>
            </div>
          </div>

          {/* Logs Panel */}
          <div className="bg-slate-900 rounded-xl shadow-inner border border-slate-800 p-4 font-mono text-xs text-green-400 h-64 flex flex-col">
            <div className="flex justify-between items-center mb-2 border-b border-slate-700 pb-2">
              <span className="flex items-center text-slate-400"><Activity className="w-4 h-4 mr-2" /> Pipeline Status</span>
              <span className="text-slate-500 text-[10px]">v2.0.0-CYP</span>
            </div>
            <div className="flex-1 overflow-y-auto space-y-1 pr-2 custom-scrollbar">
              {state.logs.map((log, i) => (
                <div key={i} className="break-words">{">"} {log}</div>
              ))}
              {state.error && <div className="text-red-400">{">"} ERROR: {state.error}</div>}
              <div ref={logsEndRef} />
            </div>
          </div>
        </div>

        {/* Right Panel: Results */}
        <div className="lg:col-span-8 flex flex-col h-[calc(100vh-8rem)] overflow-y-auto space-y-6 pr-2 custom-scrollbar">
          
          {/* Overview Card */}
          <div className="bg-white dark:bg-slate-800 rounded-xl shadow-sm border border-slate-200 dark:border-slate-700 p-6">
            <h2 className="text-xl font-bold mb-4 flex items-center text-slate-800 dark:text-slate-100"><Info className="w-5 h-5 mr-2 text-blue-500" /> Overview</h2>
            <div className="space-y-6">
              {geneInfo ? (
                    <div className="bg-blue-50 dark:bg-blue-900/20 rounded-lg p-5 border border-blue-100 dark:border-blue-800">
                      <h3 className="text-lg font-bold text-blue-900 dark:text-blue-100 mb-2">{geneInfo.symbol}</h3>
                      <p className="text-blue-800 dark:text-blue-200">{geneInfo.name}</p>
                      <div className="mt-3 flex space-x-4 text-sm">
                        <span className="bg-white dark:bg-slate-800 px-2 py-1 rounded shadow-sm text-slate-600 dark:text-slate-300">UniProt: {geneInfo.uniprot_id || 'N/A'}</span>
                        <span className="bg-white dark:bg-slate-800 px-2 py-1 rounded shadow-sm text-slate-600 dark:text-slate-300">Entrez: {geneInfo.entrez_id || 'N/A'}</span>
                      </div>
                    </div>
                  ) : (
                    <div className="text-center text-slate-500 py-10">Run analysis to see gene overview.</div>
                  )}

                  {variant && (
                    <div className="bg-white dark:bg-slate-800 rounded-lg p-5 border border-slate-200 dark:border-slate-700 shadow-sm">
                      <h3 className="text-lg font-bold mb-4 flex items-center"><Activity className="w-5 h-5 mr-2 text-purple-500" /> Variant Details: {variant.proteinChange}</h3>
                      <div className="grid grid-cols-2 gap-4">
                        <div className="p-3 bg-slate-50 dark:bg-slate-900 rounded-lg">
                          <span className="block text-xs text-slate-500 uppercase tracking-wider mb-1">Clinical Significance</span>
                          <span className={`font-medium ${variant.clinicalSignificance === 'Pathogenic' ? 'text-red-500' : variant.clinicalSignificance === 'Benign' ? 'text-green-500' : 'text-amber-500'}`}>
                            {variant.clinicalSignificance || 'Unknown'}
                          </span>
                        </div>
                        <div className="p-3 bg-slate-50 dark:bg-slate-900 rounded-lg">
                          <span className="block text-xs text-slate-500 uppercase tracking-wider mb-1">AlphaMissense Score</span>
                          <span className="font-medium">{variant.amScore !== null && variant.amScore !== undefined && !isNaN(Number(variant.amScore)) ? Number(variant.amScore).toFixed(2) : 'N/A'}</span>
                        </div>
                        <div className="p-3 bg-slate-50 dark:bg-slate-900 rounded-lg">
                          <span className="block text-xs text-slate-500 uppercase tracking-wider mb-1">ClinVar ID</span>
                          <span className="font-medium">{variant.clinVarId || 'Not Found'}</span>
                        </div>
                        <div className="p-3 bg-slate-50 dark:bg-slate-900 rounded-lg flex justify-between items-center">
                          <div>
                            <span className="block text-xs text-slate-500 uppercase tracking-wider mb-1">Target Drug</span>
                            <span className="font-medium">{drugName || 'N/A'}</span>
                          </div>
                          {drugName && geneInfo?.uniprot_id && (
                            <button 
                              onClick={handleCalculateDistances}
                              className="text-xs bg-blue-100 dark:bg-blue-900/40 text-blue-700 dark:text-blue-300 px-2 py-1 rounded hover:bg-blue-200 dark:hover:bg-blue-800/60 transition-colors"
                              title="Calculate distance from variant to drug"
                            >
                              Calc Dist
                            </button>
                          )}
                        </div>
                        <div className="p-3 bg-slate-50 dark:bg-slate-900 rounded-lg">
                          <span className="block text-xs text-slate-500 uppercase tracking-wider mb-1">GnomAD Frequency</span>
                          <span className="font-medium">{variant.gnomadFreq !== null && variant.gnomadFreq !== undefined && !isNaN(Number(variant.gnomadFreq)) ? Number(variant.gnomadFreq).toExponential(2) : 'N/A'}</span>
                        </div>
                        {variant.distanceToDrug !== undefined && variant.distanceToDrug !== null && (
                          <div className="p-3 bg-slate-50 dark:bg-slate-900 rounded-lg">
                            <span className="block text-xs text-slate-500 uppercase tracking-wider mb-1">Distance to Drug</span>
                            <span className="font-medium">{variant.distanceToDrug.toFixed(1)} Å</span>
                          </div>
                        )}
                      </div>
                    </div>
                  )}
                </div>
              </div>

              {/* Structure Card */}
              <div className="bg-white dark:bg-slate-800 rounded-xl shadow-sm border border-slate-200 dark:border-slate-700 p-6">
                <h2 className="text-xl font-bold mb-4 flex items-center text-slate-800 dark:text-slate-100"><Box className="w-5 h-5 mr-2 text-indigo-500" /> 3D Structure</h2>
                <div className="h-[600px] flex flex-col">
                  {geneInfo?.uniprot_id ? (
                    <div className="flex-1 border border-slate-200 dark:border-slate-700 rounded-lg overflow-hidden relative bg-black">
                      <StructureViewer 
                        ref={structureViewerRef}
                        humanUniprot={geneInfo.uniprot_id}
                        drugName={drugName}
                        customHighlights={variant ? [{ residue: variant.residue, color: 'red' }] : []}
                      />
                    </div>
                  ) : (
                    <div className="flex-1 flex items-center justify-center text-slate-500 border border-dashed border-slate-300 dark:border-slate-700 rounded-lg">
                      Run analysis to load 3D structure.
                    </div>
                  )}
                </div>
              </div>

              {/* Plan Card */}
              <div className="bg-white dark:bg-slate-800 rounded-xl shadow-sm border border-slate-200 dark:border-slate-700 p-6">
                <h2 className="text-xl font-bold mb-4 flex items-center text-slate-800 dark:text-slate-100"><FileText className="w-5 h-5 mr-2 text-purple-500" /> Experimental Plan</h2>
                <div className="prose prose-slate dark:prose-invert max-w-none">
                  {aiPlan ? (
                    <ReactMarkdown>{aiPlan}</ReactMarkdown>
                  ) : (
                    <div className="text-center text-slate-500 py-10 flex flex-col items-center">
                      <p className="mb-4">Experimental plan has not been generated yet.</p>
                      <button 
                        onClick={handleGeneratePlan}
                        disabled={isGeneratingPlan || !variant || !geneInfo}
                        className="bg-purple-600 hover:bg-purple-700 text-white font-medium py-2 px-6 rounded-lg transition-colors flex items-center disabled:opacity-50 shadow-sm"
                      >
                        {isGeneratingPlan ? (
                          <><Loader2 className="w-5 h-5 mr-2 animate-spin" /> Generating Plan...</>
                        ) : (
                          <><Sparkles className="w-5 h-5 mr-2" /> Generate Experimental Plan</>
                        )}
                      </button>
                    </div>
                  )}
                </div>
              </div>

              {/* CRISPR Card */}
              <div className="bg-white dark:bg-slate-800 rounded-xl shadow-sm border border-slate-200 dark:border-slate-700 p-6">
                <h2 className="text-xl font-bold mb-4 flex items-center text-slate-800 dark:text-slate-100"><Dna className="w-5 h-5 mr-2 text-emerald-500" /> CRISPR Design</h2>
                <div className="space-y-6">
                  {crisprResults.length > 0 ? (
                    <>
                      <div className="flex justify-between items-center bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                        <span className="font-medium">Found {crisprResults.length} potential CRISPR strategies</span>
                        <div className="flex space-x-2">
                          <button 
                            onClick={() => setSelectedCrisprIndex(Math.max(0, selectedCrisprIndex - 1))}
                            disabled={selectedCrisprIndex === 0}
                            className="p-1.5 bg-white dark:bg-slate-700 rounded shadow-sm disabled:opacity-50"
                          ><ChevronLeft className="w-4 h-4" /></button>
                          <span className="px-3 py-1.5 text-sm font-medium">Option {selectedCrisprIndex + 1} of {crisprResults.length}</span>
                          <button 
                            onClick={() => setSelectedCrisprIndex(Math.min(crisprResults.length - 1, selectedCrisprIndex + 1))}
                            disabled={selectedCrisprIndex === crisprResults.length - 1}
                            className="p-1.5 bg-white dark:bg-slate-700 rounded shadow-sm disabled:opacity-50"
                          ><ChevronRight className="w-4 h-4" /></button>
                        </div>
                      </div>

                      {crisprResults[selectedCrisprIndex] && (
                        <div className="space-y-6">
                          {/* Strategy Overview */}
                          <div className="grid grid-cols-2 gap-4">
                            <div className="bg-white dark:bg-slate-800 p-4 rounded-lg border border-slate-200 dark:border-slate-700 shadow-sm">
                              <span className="block text-xs text-slate-500 uppercase mb-1">Strategy</span>
                              <span className="font-medium text-blue-600 dark:text-blue-400">{crisprResults[selectedCrisprIndex].strategy.replace(/_/g, ' ')}</span>
                            </div>
                            <div className="bg-white dark:bg-slate-800 p-4 rounded-lg border border-slate-200 dark:border-slate-700 shadow-sm">
                              <span className="block text-xs text-slate-500 uppercase mb-1">Silent Mutations</span>
                              <span className="font-medium">{crisprResults[selectedCrisprIndex].silentMutationCount} introduced</span>
                            </div>
                          </div>

                          {/* Oligos */}
                          <div className="bg-white dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 shadow-sm overflow-hidden">
                            <div className="bg-slate-50 dark:bg-slate-800/80 px-4 py-3 border-b border-slate-200 dark:border-slate-700 flex justify-between items-center">
                              <h4 className="font-semibold flex items-center"><Dna className="w-4 h-4 mr-2 text-purple-500" /> Synthetic Oligonucleotides</h4>
                            </div>
                            <div className="p-4 space-y-4 font-mono text-sm">
                              <div>
                                <span className="text-xs text-slate-500 uppercase block mb-1">sgRNA Top Strand (Cloning)</span>
                                <div className="bg-slate-100 dark:bg-slate-900 p-2 rounded break-all">{crisprResults[selectedCrisprIndex].cloningOligoA}</div>
                              </div>
                              <div>
                                <span className="text-xs text-slate-500 uppercase block mb-1">sgRNA Bottom Strand (Cloning)</span>
                                <div className="bg-slate-100 dark:bg-slate-900 p-2 rounded break-all">{crisprResults[selectedCrisprIndex].cloningOligoB}</div>
                              </div>
                              <div>
                                <span className="text-xs text-slate-500 uppercase block mb-1">HDR Repair Template (ssODN)</span>
                                <div className="bg-slate-100 dark:bg-slate-900 p-2 rounded break-all text-blue-600 dark:text-blue-400">{crisprResults[selectedCrisprIndex].repairTemplate}</div>
                              </div>
                            </div>
                          </div>

                          {/* Alignment */}
                          <div className="bg-white dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 shadow-sm overflow-hidden">
                            <div className="bg-slate-50 dark:bg-slate-800/80 px-4 py-3 border-b border-slate-200 dark:border-slate-700">
                              <h4 className="font-semibold flex items-center"><List className="w-4 h-4 mr-2 text-emerald-500" /> Sequence Alignment</h4>
                            </div>
                            <div className="p-4 overflow-x-auto font-mono text-xs whitespace-pre">
                              <div className="flex"><span className="w-12 text-slate-400">REF:</span>{renderRefDna(crisprResults[selectedCrisprIndex].dnaAlignment.original, crisprResults[selectedCrisprIndex].site, crisprResults[selectedCrisprIndex].homologyStart)}</div>
                              <div className="flex"><span className="w-12 text-slate-400">VAR:</span>{renderAltSeq(crisprResults[selectedCrisprIndex].dnaAlignment.original, crisprResults[selectedCrisprIndex].dnaAlignment.modified)}</div>
                            </div>
                          </div>
                        </div>
                      )}
                    </>
                  ) : (
                    <div className="text-center text-slate-500 py-10 flex flex-col items-center">
                      <p className="mb-4">CRISPR designs have not been generated yet.</p>
                      <button 
                        onClick={handleGenerateCrispr}
                        disabled={isGeneratingCrispr || !variant || !cdna}
                        className="bg-blue-600 hover:bg-blue-700 text-white font-medium py-2 px-6 rounded-lg transition-colors flex items-center disabled:opacity-50 shadow-sm"
                      >
                        {isGeneratingCrispr ? (
                          <><Loader2 className="w-5 h-5 mr-2 animate-spin" /> Generating Designs...</>
                        ) : (
                          <><Dna className="w-5 h-5 mr-2" /> Generate CRISPR Designs</>
                        )}
                      </button>
                    </div>
                  )}
                </div>
              </div>

              {/* Modelling Card */}
              <div className="bg-white dark:bg-slate-800 rounded-xl shadow-sm border border-slate-200 dark:border-slate-700 p-6">
                <h2 className="text-xl font-bold mb-4 flex items-center text-slate-800 dark:text-slate-100"><Activity className="w-5 h-5 mr-2 text-pink-500" /> Molecular Modelling (DiffDock)</h2>
                <div className="space-y-6">
                  {diffDockResults ? (
                    <div className="bg-slate-50 dark:bg-slate-900 p-4 rounded-lg border border-slate-200 dark:border-slate-700 overflow-x-auto">
                      <pre className="text-xs font-mono text-slate-800 dark:text-slate-300">
                        {JSON.stringify(diffDockResults, null, 2)}
                      </pre>
                    </div>
                  ) : (
                    <div className="text-center text-slate-500 py-10 flex flex-col items-center">
                      <p className="mb-4">Molecular modelling has not been run yet.</p>
                      <button 
                        onClick={handleGenerateDiffDock}
                        disabled={isGeneratingDiffDock}
                        className="bg-pink-600 hover:bg-pink-700 text-white font-medium py-2 px-6 rounded-lg transition-colors flex items-center disabled:opacity-50 shadow-sm"
                      >
                        {isGeneratingDiffDock ? (
                          <><Loader2 className="w-5 h-5 mr-2 animate-spin" /> Running DiffDock...</>
                        ) : (
                          <><Activity className="w-5 h-5 mr-2" /> Run DiffDock</>
                        )}
                      </button>
                    </div>
                  )}
                </div>
              </div>

        </div>
      </main>
    </div>
  );
};

export default App;
