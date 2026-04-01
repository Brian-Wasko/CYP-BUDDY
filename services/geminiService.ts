
import { GoogleGenAI } from "@google/genai";
import { Variant, TargetCYP } from '../types';

export const generateExperimentalPlan = async (
  humanGene: TargetCYP,
  drugName: string,
  variant: Variant,
  onUpdate: (text: string) => void,
  structureImage?: string | null
): Promise<void> => {
  
  const variantContext = `Variant: ${variant.proteinChange} (Ref: ${variant.refAA} -> Mut: ${variant.targetAA}). AlphaMissense Score: ${variant.amScore ?? 'N/A'}. Clinical Significance: ${variant.clinicalSignificance || 'Unknown'}.`;

  let imageContext = "";
  if (structureImage) {
      imageContext = `\n\nATTACHED IMAGE: A snapshot of the 3D protein structure showing the variant location.`;
  }

  const prompt = `
    Gene: ${humanGene} (Human Cytochrome P450)
    Drug: ${drugName}
    ${variantContext}
    ${imageContext}

    Task: Act as a grounded research assistant. Communicate like a scientist. Use Web Search (grounding) to find relevant papers regarding the human CYP gene, the specific variant, and the drug.

    **STEP 2: Literature Search & Drug Interaction**
    - Identify known interactions between ${humanGene} and the drug "${drugName}".
    - Explicitly state if peer-reviewed literature suggests this specific SNP (${variant.proteinChange}) alters the metabolism of the drug. 
    - If literature is absent for this specific SNP-drug interaction, state clearly that this is a novel discovery opportunity for the student.
    - **Crucially, suggest expected CYP-dependent growth phenotypes in the yeast assay based on the known metabolism of ${drugName} by ${humanGene} (e.g., does the CYP bioactivate a prodrug into a toxic metabolite, leading to increased sensitivity, or does it detoxify a drug, leading to resistance?).**

    **STEP 3: Experimental Design (The 12Δ "Intracellular Trap" Assay)**
    - Design a 96-well or 384-well IC50 growth inhibition assay.
    - **Control Group:** The parental 12Δ hypersensitive yeast strain (lacking the human CYP and lacking 12 key drug efflux pumps).
    - **Experimental Group:** The engineered 12Δ strain expressing the wild-type ${humanGene}, alongside the student's CRISPR-edited 12Δ strain expressing their specific ${humanGene} variant (${variant.proteinChange}).
    - Emphasize that monitoring shifted IC50 values (bioactivation to a toxic metabolite, or detoxification) replaces the need for animal testing in early-stage toxicology.

    **CITATION REQUIREMENT:** Include complete citations whenever you reference specific findings or protocols that were found via web search.

    Format Rules:
    - Use clean Markdown.
    - SECTION TITLES MUST BE IN ALL CAPS.
    - Insert a horizontal rule (---) between each section.
    - Always end with a last sentence: "Disclaimer: This content is AI-generated. Strongly consider validating phenotype(s) and other provided information using primary literature. All experiments should be reviewed and overseen by a qualified scientist for safety and compliance with environmental health and safety regulations."
  `;

  try {
    const ai = new GoogleGenAI({ apiKey: process.env.GEMINI_API_KEY });
    
    const parts: any[] = [{ text: prompt }];
    if (structureImage) {
        const match = structureImage.match(/^data:(image\/[a-zA-Z]+);base64,(.+)$/);
        if (match) {
            parts.push({
                inlineData: {
                    mimeType: match[1],
                    data: match[2]
                }
            });
        } else {
            const base64Data = structureImage.split(',')[1] || structureImage;
            parts.push({
                inlineData: {
                    mimeType: "image/png",
                    data: base64Data
                }
            });
        }
    }

    const response = await ai.models.generateContent({
      model: 'gemini-3.1-pro-preview', 
      contents: { parts },
      config: {
        tools: [{ googleSearch: {} }],
      },
    });

    let fullText = response.text || "";
    let collectedGroundingChunks: any[] = [];

    if (response.candidates?.[0]?.groundingMetadata?.groundingChunks) {
      collectedGroundingChunks = response.candidates[0].groundingMetadata.groundingChunks;
    }

    if (collectedGroundingChunks.length > 0) {
      const validLinks = collectedGroundingChunks
        .filter((c: any) => c.web?.uri && c.web?.title)
        .map((c: any) => `- [${c.web.title}](${c.web.uri})`);
      
      if (validLinks.length > 0) {
        fullText += `\n\n---\n### REFERENCES & SOURCES\n${validLinks.join('\n')}`;
      }
    }

    onUpdate(fullText);

  } catch (error: any) {
    console.error("Gemini Error:", error);
    if (error.message?.includes("429") || error.status === 429 || error.toString().includes("RESOURCE_EXHAUSTED")) {
         throw new Error("Gemini API Quota Exceeded (429). Please wait a minute and try again, or use a paid API key.");
    }
    throw error;
  }
};

export const fetchSuggestedDrugs = async (humanGene: TargetCYP): Promise<string[]> => {
  try {
    const ai = new GoogleGenAI({ apiKey: process.env.GEMINI_API_KEY });
    const response = await ai.models.generateContent({
      model: 'gemini-3-flash-preview',
      contents: `List 10 FDA-approved drugs that are primarily metabolized by human ${humanGene}. Return ONLY a JSON array of strings.`,
      config: {
        responseMimeType: "application/json",
      },
    });
    
    const text = response.text || "[]";
    return JSON.parse(text);
  } catch (error) {
    console.error("Error fetching suggested drugs:", error);
    return [];
  }
};

