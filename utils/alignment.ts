
export const AA_MAP: Record<string, string> = {
  'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
  'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
  'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
  'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
  'Ter': '*'
};

export const parseProteinChange = (pChange: string): { ref: string, res: number, target: string } | null => {
  // Format p.Arg114Gln
  const regex3 = /p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})/;
  const match3 = pChange.match(regex3);
  if (match3) {
    return {
      ref: AA_MAP[match3[1]] || '?',
      res: parseInt(match3[2]),
      target: AA_MAP[match3[3]] || '?'
    };
  }

  // Format R144C
  const regex1 = /^([A-Z])(\d+)([A-Z\*])$/i;
  const match1 = pChange.match(regex1);
  if (match1) {
    return {
      ref: match1[1].toUpperCase(),
      res: parseInt(match1[2]),
      target: match1[3].toUpperCase()
    };
  }

  return null;
};
