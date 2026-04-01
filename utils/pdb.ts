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

export const fetchPdbForGeneAndDrug = async (uniprotId: string, drugName: string): Promise<string | null> => {
    const cleanId = uniprotId.trim().toUpperCase();
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
    
    try {
        const rcsbRes = await fetch(`https://search.rcsb.org/rcsbsearch/v2/query?json=${encodeURIComponent(JSON.stringify(combinedQuery))}`);
        if (rcsbRes.ok) {
            const rcsbData = await rcsbRes.json();
            if (rcsbData.result_set && rcsbData.result_set.length > 0) {
                const pdbId = rcsbData.result_set[0].identifier;
                return await fetchWithProxy(`https://files.rcsb.org/download/${pdbId}.pdb`);
            }
        }
    } catch (e) {
        console.warn("Failed to fetch PDB for distance calculation", e);
    }
    return null;
};

export const calculateDistancesToDrug = (pdbData: string, residues: number[]): Record<number, number> => {
    const atoms: { resi: number, x: number, y: number, z: number, isHet: boolean, resn: string }[] = [];
    const lines = pdbData.split('\n');
    for (const line of lines) {
        if (line.startsWith('ATOM  ') || line.startsWith('HETATM')) {
            const isHet = line.startsWith('HETATM');
            const resn = line.substring(17, 20).trim();
            if (isHet && resn === 'HOH') continue; // Skip water
            
            const resi = parseInt(line.substring(22, 26).trim(), 10);
            const x = parseFloat(line.substring(30, 38).trim());
            const y = parseFloat(line.substring(38, 46).trim());
            const z = parseFloat(line.substring(46, 54).trim());
            
            if (!isNaN(resi) && !isNaN(x) && !isNaN(y) && !isNaN(z)) {
                atoms.push({ resi, x, y, z, isHet, resn });
            }
        }
    }

    const drugAtoms = atoms.filter(a => a.isHet);
    const proteinAtoms = atoms.filter(a => !a.isHet);

    const distances: Record<number, number> = {};

    if (drugAtoms.length === 0) return distances;

    for (const resi of residues) {
        const resAtoms = proteinAtoms.filter(a => a.resi === resi);
        if (resAtoms.length === 0) continue;

        let minDistance = Infinity;
        for (const a1 of resAtoms) {
            for (const a2 of drugAtoms) {
                const dx = a1.x - a2.x;
                const dy = a1.y - a2.y;
                const dz = a1.z - a2.z;
                const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                if (dist < minDistance) minDistance = dist;
            }
        }
        distances[resi] = minDistance;
    }

    return distances;
};
