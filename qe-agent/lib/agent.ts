import { ChatOpenAI } from "@langchain/openai";
import { createReactAgent } from "@langchain/langgraph/prebuilt";
import { 
  generateQEInputTool, 
  validateQEParametersTool, 
  suggestCalculationTool 
} from "./qe-tools";

// Initialize the LLM
export function createQEAgent() {
  const model = new ChatOpenAI({
    modelName: process.env.OPENAI_MODEL || "gpt-4o-mini",
    temperature: 0.1,
    openAIApiKey: process.env.OPENAI_API_KEY,
  });

  // Create agent with QE tools
  const tools = [
    generateQEInputTool,
    validateQEParametersTool,
    suggestCalculationTool,
  ];

  const agent = createReactAgent({
    llm: model,
    tools,
  });

  return agent;
}

// System message for the agent
export const QE_SYSTEM_MESSAGE = `You are an expert assistant for Quantum ESPRESSO (QE), a software package for electronic structure calculations.

Your role is to help users create properly formatted QE input files. You have access to tools that can:
1. Generate QE input files with proper formatting
2. Validate parameters and provide recommendations
3. Suggest appropriate calculation types

When a user asks to create a QE input file:
1. Ask clarifying questions if needed (calculation type, material, parameters)
2. Use the suggest_calculation_type tool if the user is unsure about calculation type
3. Validate parameters using the validate_qe_parameters tool
4. Generate the input file using the generate_qe_input tool
5. Provide the complete input file to the user

Common calculation types:
- scf: Self-consistent field calculation (ground state)
- relax: Geometry optimization (fixed cell)
- vc-relax: Variable-cell relaxation (optimize cell and positions)
- nscf: Non-self-consistent calculation (for DOS, bands)
- bands: Band structure calculation

Important parameters:
- ecutwfc: Wavefunction cutoff energy (Ry) - typically 30-80 Ry
- ecutrho: Charge density cutoff (Ry) - typically 4-12 times ecutwfc
- k_points: K-point sampling for Brillouin zone integration

Be helpful, ask for missing information, and provide explanations when needed.`;
