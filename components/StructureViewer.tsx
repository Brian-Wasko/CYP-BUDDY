import React, { useEffect, useRef, useState, useMemo } from 'react';
import { Camera, ExternalLink, AlertCircle, Box, Palette, Database, Save, X, Dna } from 'lucide-react';

declare global {
  interface Window {
    $3Dmol: any;
  }
}

interface Highlight {
  residue: number;
  color: string;
}

interface Props {
  humanUniprot: string | null;
  customHighlights: Highlight[];
  drugName?: string;
}

type Representation = 'cartoon' | 'stick' | 'surface' | 'sphere';

export interface StructureViewerHandle {
  captureImage: () => string | null;
}

export const StructureViewer = React.forwardRef<StructureViewerHandle, Props>(({ 
  humanUniprot, 
  customHighlights,
  drugName
}, ref) => {
  const [representation, setRepresentation] = useState<Representation>('cartoon');
  const [variantRepresentation, setVariantRepresentation] = useState<Representation>('sphere');
  
  // ID State (allows manual override)
  const [activeHumanId, setActiveHumanId] = useState<string | null>(humanUniprot);
  const [showIdControls, setShowIdControls] = useState(false);
  const [humanInput, setHumanInput] = useState('');

  // Expose captureImage method
  React.useImperativeHandle(ref, () => ({
    captureImage: () => {
        if (!viewerRef.current) return null;
        const v = viewerRef.current;
        v.render();
        return v.pngURI();
    }
  }));

  // Colors State
  const [colors, setColors] = useState({
      singleBase: '#94a3b8', // Slate-400
      singleVariant: '#ef4444', // Red-500
      drug: '#10b981', // Emerald-500
      heme: '#ea580c', // Orange-500
  });

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [structureInfo, setStructureInfo] = useState<{ source: string, id: string } | null>(null);
  const [proteinFeatures, setProteinFeatures] = useState<any[]>([]);
  
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);

  // Sync props to state when gene changes
  useEffect(() => {
      setActiveHumanId(humanUniprot);
      setHumanInput(humanUniprot || '');
      if (humanUniprot) {
          fetchProteinFeatures(humanUniprot);
      } else {
          setProteinFeatures([]);
      }
  }, [humanUniprot]);

  const fetchProteinFeatures = async (uniprotId: string) => {
    try {
      const res = await fetch(`https://rest.uniprot.org/uniprotkb/${uniprotId}`);
      if (res.ok) {
        const data = await res.json();
        if (data.features) {
          setProteinFeatures(data.features);
        }
      }
    } catch (e) {
      console.warn("Failed to fetch UniProt features", e);
    }
  };

  const fetchWithProxy = async (url: string) => {
    try {
        const res = await fetch(url);
        if (res.ok) return await res.text();
    } catch (e) {
        // Fallback
    }
    try {
        const res = await fetch(`https://corsproxy.io/?${encodeURIComponent(url)}`);
        if (res.ok) return await res.text();
    } catch (e) {
        console.warn("CorsProxy failed", e);
    }
    try {
        const res = await fetch(`https://api.allorigins.win/get?url=${encodeURIComponent(url)}`);
        if (res.ok) {
            const json = await res.json();
            return json.contents;
        }
    } catch (e) {
        console.warn("AllOrigins failed", e);
    }
    throw new Error("Failed to fetch structure data.");
  };

  const getPdbData = async (id: string): Promise<{ data: string, source: string, id: string }> => {
      const cleanId = id.trim().toUpperCase();
      
      // Heuristic: PDB IDs are exactly 4 chars. UniProt IDs are usually 6+.
      const isLikelyPdb = cleanId.length === 4;

      if (isLikelyPdb) {
          try {
              const data = await fetchWithProxy(`https://files.rcsb.org/download/${cleanId}.pdb`);
              return { data, source: 'PDB(RCSB)', id: cleanId };
          } catch (e) {
              console.warn("RCSB direct fetch failed, trying AlphaFold as fallback");
          }
      }

      // If we have a drug name, try to find a co-complex in PDB first
      if (drugName && drugName.trim().length > 0) {
          try {
              const combinedQuery = {
                  query: {
                      type: "group",
                      logical_operator: "and",
                      nodes: [
                          {
                              type: "terminal",
                              service: "text",
                              parameters: {
                                  attribute: "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                                  operator: "exact_match",
                                  value: cleanId
                              }
                          },
                          {
                              type: "terminal",
                              service: "full_text",
                              parameters: {
                                  value: drugName.trim()
                              }
                          }
                      ]
                  },
                  return_type: "entry",
                  request_options: {
                      scoring_strategy: "combined",
                      sort: [{ sort_by: "score", direction: "desc" }]
                  }
              };
              
              const rcsbRes = await fetch(`https://search.rcsb.org/rcsbsearch/v2/query?json=${encodeURIComponent(JSON.stringify(combinedQuery))}`);
              if (rcsbRes.ok) {
                  const rcsbData = await rcsbRes.json();
                  if (rcsbData.result_set && rcsbData.result_set.length > 0) {
                      const pdbId = rcsbData.result_set[0].identifier;
                      const data = await fetchWithProxy(`https://files.rcsb.org/download/${pdbId}.pdb`);
                      return { data, source: `RCSB (with ${drugName})`, id: pdbId };
                  }
              }
          } catch (e) {
              console.warn("RCSB combined fetch failed", e);
          }
      }

      // Try AlphaFold first for UniProt IDs if no drug co-complex found
      try {
        const afRes = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${cleanId}`);
        if (afRes.ok) {
            const afData = await afRes.json();
            if (afData && afData.length > 0) {
                const pdbUrl = afData[0].pdbUrl;
                if (pdbUrl) {
                    const data = await fetchWithProxy(pdbUrl);
                    return { data, source: 'AlphaFold', id: cleanId };
                }
            }
        }
      } catch (e) {
          console.warn("AlphaFold fetch failed", e);
      }

      // Fallback to RCSB search if AlphaFold failed or if it wasn't a PDB ID initially
      try {
          const query = {
              query: {
                  type: "terminal",
                  service: "text",
                  parameters: {
                      attribute: "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                      operator: "exact_match",
                      value: cleanId
                  }
              },
              return_type: "entry",
              request_options: {
                  scoring_strategy: "combined",
                  sort: [{ sort_by: "score", direction: "desc" }]
              }
          };
          
          const rcsbRes = await fetch(`https://search.rcsb.org/rcsbsearch/v2/query?json=${encodeURIComponent(JSON.stringify(query))}`);
          if (rcsbRes.ok) {
              const rcsbData = await rcsbRes.json();
              if (rcsbData.result_set && rcsbData.result_set.length > 0) {
                  const pdbId = rcsbData.result_set[0].identifier;
                  const data = await fetchWithProxy(`https://files.rcsb.org/download/${pdbId}.pdb`);
                  return { data, source: 'RCSB (Lookup)', id: pdbId };
              }
          }
      } catch (e) {
          console.warn("RCSB fetch failed", e);
      }

      throw new Error(`No structure found for ${cleanId}`);
  };

  const getStyle = (rep: Representation, color: string, opacity: number = 1.0) => {
      if (rep === 'cartoon') return { cartoon: { color, opacity } };
      if (rep === 'stick') return { stick: { color, opacity, radius: 0.2 } };
      if (rep === 'sphere') return { sphere: { color, opacity, scale: 1.0 } };
      // For surface, we usually use the model as base or hide it
      return { cartoon: { color, opacity: 0 } }; 
  };

  const loadStructure = async () => {
    setLoading(true);
    setError(null);
    setStructureInfo(null);

    if (!viewerRef.current && containerRef.current && window.$3Dmol) {
        const config = { backgroundColor: '#1e293b' }; // slate-800
        viewerRef.current = window.$3Dmol.createViewer(containerRef.current, config);
    }
    
    const v = viewerRef.current;
    if (!v) {
        setLoading(false);
        return;
    }
    
    v.clear();

    try {
        if (!activeHumanId) throw new Error(`No ID provided`);
        
        const res = await getPdbData(activeHumanId);
        const model = v.addModel(res.data, "pdb");
        
        v.custom_models = { single: model };
        
        setStructureInfo({ source: res.source, id: res.id });

        v.zoomTo();
        updateRender(); // Apply initial styles

    } catch (err) {
        setError((err as Error).message);
    } finally {
        setLoading(false);
    }
  };

  // Function to apply styles based on current state (colors, representation) without reloading geometry
  const updateRender = () => {
      const v = viewerRef.current;
      if (!v || !v.custom_models) return;

      v.removeAllSurfaces();
      v.removeAllLabels();
      v.setStyle({}, {}); // Clear all existing styles

      const applyVariantStyle = (model: any, residue: number, color: string) => {
          const sel = { model: model, resi: residue };
          
          if (variantRepresentation === 'surface') {
              v.addSurface(window.$3Dmol.SurfaceType.VDW, { opacity: 1.0, color: color }, sel);
          } else {
              let style: any = {};
              if (variantRepresentation === 'stick') style = { stick: { color, radius: 0.3 } };
              else if (variantRepresentation === 'sphere') style = { sphere: { color, scale: 0.8 } };
              else if (variantRepresentation === 'cartoon') style = { cartoon: { color, thickness: 1.0, opacity: 1.0 } };
              
              v.addStyle(sel, style);
          }
      };

      // Single Mode
      const model = v.custom_models.single;
      
      // Safety check
      if (!model) return;

      const highlights = customHighlights || [];

      if (representation === 'surface') {
          v.addSurface(window.$3Dmol.SurfaceType.VDW, {opacity: 0.9, color: colors.singleBase}, {model: model, hetflag: false});
          model.setStyle({hetflag: true}, { stick: { color: colors.drug, radius: 0.3 } });
      } else {
          model.setStyle({hetflag: false}, getStyle(representation, colors.singleBase));
          model.setStyle({hetflag: true}, { stick: { color: colors.drug, radius: 0.3 } });
      }
      model.setStyle({resn: 'HEM'}, { stick: { color: colors.heme, radius: 0.3 } }); // Specific style for Heme
      model.setStyle({resn: 'HOH'}, { hidden: true }); // hide water

      // Apply protein features (Active sites & Binding sites)
      proteinFeatures.forEach(feature => {
          if (feature.type === 'Active site' || feature.type === 'Binding site') {
              const isSubstrate = feature.description?.toLowerCase().includes('substrate') || feature.type === 'Binding site';
              const color = feature.type === 'Active site' ? '#3b82f6' : (isSubstrate ? '#eab308' : null); // blue for active, yellow for binding
              
              if (color) {
                  let start = feature.location?.start?.value || feature.location?.position?.value;
                  let end = feature.location?.end?.value || feature.location?.position?.value;
                  
                  if (start && end) {
                      for (let i = start; i <= end; i++) {
                          const sel = { model: model, resi: i };
                          // Highlight the cartoon and add sticks
                          if (representation === 'cartoon') {
                              model.setStyle(sel, { cartoon: { color: color } });
                              v.addStyle(sel, { stick: { color: color, radius: 0.15 } });
                          } else if (representation === 'surface') {
                              v.addSurface(window.$3Dmol.SurfaceType.VDW, {opacity: 1.0, color: color}, sel);
                          } else {
                              model.setStyle(sel, getStyle(representation, color));
                          }
                      }
                  }
              }
          }
      });

      // Variants
      highlights?.forEach((h: Highlight) => {
          applyVariantStyle(model, h.residue, colors.singleVariant);
          // Label
          v.addLabel(`${h.residue}`, { 
              position: { resi: h.residue }, 
              backgroundColor: 'rgba(0,0,0,0.7)', 
              fontColor: 'white', 
              fontSize: 12 
          });
      });
      
      v.render();
  };

  useEffect(() => {
      if (window.$3Dmol) {
          loadStructure();
      } else {
          setError("3Dmol.js library not loaded.");
      }
  }, [activeHumanId, drugName]); 

  // Re-render when display options change
  useEffect(() => {
      updateRender();
  }, [representation, variantRepresentation, colors, customHighlights, proteinFeatures]);

  const handleSnapshot = () => {
      if (!viewerRef.current) return;
      const v = viewerRef.current;
      v.render();
      const dataURI = v.pngURI();
      const link = document.createElement('a');
      link.href = dataURI;
      link.download = `Structure_${structureInfo?.id || 'protein'}.png`;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
  };

  const handleSaveIds = () => {
      setActiveHumanId(humanInput || null);
      setShowIdControls(false);
  };

  return (
    <div className="bg-slate-800 rounded-xl border border-slate-700 shadow-sm overflow-hidden flex flex-col h-[650px] animate-fade-in">
        {/* Header */}
        <div className="px-4 py-3 bg-slate-900 border-b border-slate-700 flex justify-between items-center relative">
            <div className="flex items-center gap-2">
                <Box className="w-5 h-5 text-blue-500" />
                <h3 className="font-bold text-slate-100">3D Structure</h3>
                {structureInfo && (
                    <div className="hidden sm:flex items-center gap-2">
                        <span className="text-[10px] bg-slate-700 text-slate-300 px-2 py-0.5 rounded-full font-medium flex items-center gap-1">
                            {structureInfo.source}: {structureInfo.id}
                        </span>
                        <button 
                            onClick={() => setShowIdControls(!showIdControls)}
                            className="p-1 rounded-full hover:bg-slate-700 text-slate-400 hover:text-white transition-colors"
                            title="Manually Change Structure ID (PDB/UniProt)"
                        >
                            <Database className="w-3 h-3" />
                        </button>
                    </div>
                )}
            </div>

            {/* ID Override Popup */}
            {showIdControls && (
                <div className="absolute top-12 left-4 z-50 bg-slate-800 border border-slate-600 rounded-lg p-3 shadow-xl w-64 animate-fade-in">
                    <div className="flex justify-between items-center mb-2">
                        <h4 className="text-xs font-bold text-slate-300 uppercase">Change Structure</h4>
                        <button onClick={() => setShowIdControls(false)} className="text-slate-500 hover:text-white"><X className="w-3 h-3" /></button>
                    </div>
                    <div className="space-y-2 mb-3">
                        <div>
                            <label className="block text-[10px] text-slate-500 mb-1">Human ID (UniProt or PDB)</label>
                            <input 
                                type="text" 
                                value={humanInput} 
                                onChange={(e) => setHumanInput(e.target.value)}
                                className="w-full bg-slate-900 border border-slate-700 rounded px-2 py-1 text-xs text-white focus:border-blue-500 outline-none"
                                placeholder="e.g. P12345 or 1ABC"
                            />
                        </div>
                    </div>
                    <button 
                        onClick={handleSaveIds}
                        className="w-full py-1 bg-blue-600 hover:bg-blue-700 text-white rounded text-xs font-bold flex items-center justify-center gap-1"
                    >
                        <Save className="w-3 h-3" /> Update Viewer
                    </button>
                </div>
            )}

            <div className="flex items-center gap-2">
                <button 
                    onClick={handleSnapshot}
                    disabled={loading || !!error}
                    className="p-1.5 text-slate-400 hover:text-blue-400 hover:bg-slate-700 rounded-lg transition-colors"
                    title="Take HD Snapshot"
                >
                    <Camera className="w-4 h-4" />
                </button>
            </div>
        </div>

        {/* Viewer Area */}
        <div className="flex-grow relative bg-slate-900 group">
            {loading && (
                <div className="absolute inset-0 flex flex-col items-center justify-center bg-slate-900/80 z-10 text-blue-500">
                    <Dna className="w-8 h-8 animate-spin mb-2" />
                    <span className="text-xs font-bold uppercase tracking-wider">Fetching Structure...</span>
                </div>
            )}
            
            {error && (
                <div className="absolute inset-0 flex flex-col items-center justify-center bg-slate-900 z-10 text-slate-400 p-6 text-center">
                    <AlertCircle className="w-8 h-8 mb-2 text-red-400" />
                    <p className="text-sm font-medium mb-4">{error}</p>
                    <div className="flex gap-2 text-xs justify-center">
                        {activeHumanId && (
                            <a href={`https://alphafold.ebi.ac.uk/entry/${activeHumanId}`} target="_blank" rel="noreferrer" className="text-blue-400 hover:underline flex items-center gap-1">
                                AF-Human <ExternalLink className="w-3 h-3"/>
                            </a>
                        )}
                    </div>
                    <button onClick={() => setShowIdControls(true)} className="mt-4 px-3 py-1 bg-slate-700 hover:bg-slate-600 rounded text-xs text-white">
                        Try different ID
                    </button>
                </div>
            )}

            <div id="3dmol-container" ref={containerRef} className="w-full h-full cursor-move relative z-0"></div>
        </div>

        {/* Controls Bar */}
        <div className="px-4 py-3 bg-slate-900 border-t border-slate-700 flex flex-col gap-3">
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-xs items-center">
                {/* View Style Options */}
                <div className="flex flex-col gap-2">
                    <div className="flex items-center gap-2">
                        <span className="font-bold text-slate-400 uppercase tracking-wider text-[10px] w-20">Base Style</span>
                        <div className="flex bg-slate-800 rounded border border-slate-600 p-0.5">
                            <button onClick={() => setRepresentation('cartoon')} className={`px-2 py-1 rounded ${representation === 'cartoon' ? 'bg-slate-600 text-white' : 'text-slate-400 hover:text-slate-300'}`}>Cartoon</button>
                            <button onClick={() => setRepresentation('stick')} className={`px-2 py-1 rounded ${representation === 'stick' ? 'bg-slate-600 text-white' : 'text-slate-400 hover:text-slate-300'}`}>Stick</button>
                            <button onClick={() => setRepresentation('surface')} className={`px-2 py-1 rounded ${representation === 'surface' ? 'bg-slate-600 text-white' : 'text-slate-400 hover:text-slate-300'}`}>Surface</button>
                        </div>
                    </div>
                    <div className="flex items-center gap-2">
                        <span className="font-bold text-slate-400 uppercase tracking-wider text-[10px] w-20">Variant Style</span>
                        <div className="flex bg-slate-800 rounded border border-slate-600 p-0.5">
                            <button onClick={() => setVariantRepresentation('cartoon')} className={`px-2 py-1 rounded ${variantRepresentation === 'cartoon' ? 'bg-slate-600 text-white' : 'text-slate-400 hover:text-slate-300'}`}>Cartoon</button>
                            <button onClick={() => setVariantRepresentation('stick')} className={`px-2 py-1 rounded ${variantRepresentation === 'stick' ? 'bg-slate-600 text-white' : 'text-slate-400 hover:text-slate-300'}`}>Stick</button>
                            <button onClick={() => setVariantRepresentation('sphere')} className={`px-2 py-1 rounded ${variantRepresentation === 'sphere' ? 'bg-slate-600 text-white' : 'text-slate-400 hover:text-slate-300'}`}>Sphere</button>
                            <button onClick={() => setVariantRepresentation('surface')} className={`px-2 py-1 rounded ${variantRepresentation === 'surface' ? 'bg-slate-600 text-white' : 'text-slate-400 hover:text-slate-300'}`}>Surface</button>
                        </div>
                    </div>
                </div>

                {/* Color Controls */}
                <div className="flex flex-col items-end gap-2">
                    <div className="flex items-center justify-end gap-3 flex-wrap">
                        <span className="font-bold text-slate-400 uppercase tracking-wider text-[10px] flex items-center gap-1">
                            <Palette className="w-3 h-3" /> Colors
                        </span>
                        
                        <div className="flex items-center gap-1">
                            <input type="color" value={colors.singleBase} onChange={e => setColors({...colors, singleBase: e.target.value})} className="w-6 h-6 rounded cursor-pointer bg-transparent border-0 p-0" />
                            <span className="text-slate-500">Protein</span>
                        </div>
                        <div className="flex items-center gap-1">
                            <input type="color" value={colors.singleVariant} onChange={e => setColors({...colors, singleVariant: e.target.value})} className="w-6 h-6 rounded cursor-pointer bg-transparent border-0 p-0" />
                            <span className="text-slate-500">Variant</span>
                        </div>
                        <div className="flex items-center gap-1">
                            <input type="color" value={colors.drug} onChange={e => setColors({...colors, drug: e.target.value})} className="w-6 h-6 rounded cursor-pointer bg-transparent border-0 p-0" />
                            <span className="text-slate-500">Drug</span>
                        </div>
                        <div className="flex items-center gap-1">
                            <input type="color" value={colors.heme} onChange={e => setColors({...colors, heme: e.target.value})} className="w-6 h-6 rounded cursor-pointer bg-transparent border-0 p-0" />
                            <span className="text-slate-500">Heme</span>
                        </div>
                    </div>
                    
                    {/* Feature Legend */}
                    <div className="flex items-center justify-end gap-3 text-[10px] font-medium mt-1">
                        <div className="flex items-center gap-1">
                            <div className="w-2 h-2 rounded-full bg-blue-500"></div>
                            <span className="text-slate-400">Active Site</span>
                        </div>
                        <div className="flex items-center gap-1">
                            <div className="w-2 h-2 rounded-full bg-yellow-500"></div>
                            <span className="text-slate-400">Binding Site</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
  );
});
