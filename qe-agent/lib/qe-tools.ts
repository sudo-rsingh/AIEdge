import { DynamicStructuredTool } from "@langchain/core/tools";
import { z } from "zod";

// Schema for QE input parameters
export const QEInputSchema = z.object({
  calculation: z.enum(['scf', 'relax', 'vc-relax', 'nscf', 'bands']).describe('Type of calculation'),
  prefix: z.string().describe('Prefix for output files'),
  pseudo_dir: z.string().default('./pseudo').describe('Directory containing pseudopotentials'),
  outdir: z.string().default('./out').describe('Output directory'),
  ibrav: z.number().describe('Bravais lattice index (0 for free lattice)'),
  nat: z.number().describe('Number of atoms'),
  ntyp: z.number().describe('Number of atom types'),
  ecutwfc: z.number().describe('Kinetic energy cutoff for wavefunctions (Ry)'),
  ecutrho: z.number().optional().describe('Kinetic energy cutoff for charge density (Ry)'),
  occupations: z.string().default('smearing').describe('Occupation type'),
  smearing: z.string().default('gaussian').describe('Smearing type'),
  degauss: z.number().default(0.02).describe('Smearing parameter'),
  mixing_beta: z.number().default(0.7).describe('Mixing factor for self-consistency'),
  conv_thr: z.number().default(1.0e-6).describe('Convergence threshold'),
  celldm1: z.number().optional().describe('Lattice parameter (bohr)'),
  a: z.number().optional().describe('Lattice parameter a (angstrom)'),
  b: z.number().optional().describe('Lattice parameter b (angstrom)'),
  c: z.number().optional().describe('Lattice parameter c (angstrom)'),
  atom_types: z.array(z.object({
    name: z.string(),
    mass: z.number(),
    pseudopotential: z.string()
  })).describe('Atomic species information'),
  atom_positions: z.array(z.object({
    type: z.string(),
    x: z.number(),
    y: z.number(),
    z: z.number()
  })).describe('Atomic positions'),
  k_points: z.object({
    type: z.enum(['automatic', 'gamma', 'crystal']),
    grid: z.array(z.number()).optional()
  }).describe('K-points sampling'),
});

export type QEInputParams = z.infer<typeof QEInputSchema>;

// Tool to generate QE input file
export const generateQEInputTool = new DynamicStructuredTool({
  name: "generate_qe_input",
  description: "Generate a Quantum Espresso input file based on the provided parameters. This tool creates a properly formatted PWscf input file.",
  schema: QEInputSchema,
  func: async (params: QEInputParams) => {
    return generateQEInput(params);
  },
});

// Function to generate QE input file content
export function generateQEInput(params: QEInputParams): string {
  const {
    calculation,
    prefix,
    pseudo_dir,
    outdir,
    ibrav,
    nat,
    ntyp,
    ecutwfc,
    ecutrho,
    occupations,
    smearing,
    degauss,
    mixing_beta,
    conv_thr,
    celldm1,
    a, b, c,
    atom_types,
    atom_positions,
    k_points,
  } = params;

  let input = `&CONTROL
  calculation = '${calculation}'
  prefix = '${prefix}'
  pseudo_dir = '${pseudo_dir}'
  outdir = '${outdir}'
/

&SYSTEM
  ibrav = ${ibrav}
  nat = ${nat}
  ntyp = ${ntyp}
  ecutwfc = ${ecutwfc}
`;

  if (ecutrho) {
    input += `  ecutrho = ${ecutrho}\n`;
  }

  // Add lattice parameters
  if (celldm1) {
    input += `  celldm(1) = ${celldm1}\n`;
  } else if (a) {
    input += `  A = ${a}\n`;
    if (b) input += `  B = ${b}\n`;
    if (c) input += `  C = ${c}\n`;
  }

  input += `  occupations = '${occupations}'
  smearing = '${smearing}'
  degauss = ${degauss}
/

&ELECTRONS
  mixing_beta = ${mixing_beta}
  conv_thr = ${conv_thr}
/
`;

  // Add IONS section for relaxation
  if (calculation === 'relax' || calculation === 'vc-relax') {
    input += `
&IONS
/
`;
  }

  // Add CELL section for vc-relax
  if (calculation === 'vc-relax') {
    input += `
&CELL
/
`;
  }

  // Add ATOMIC_SPECIES
  input += `\nATOMIC_SPECIES\n`;
  atom_types.forEach(type => {
    input += `  ${type.name}  ${type.mass}  ${type.pseudopotential}\n`;
  });

  // Add ATOMIC_POSITIONS
  input += `\nATOMIC_POSITIONS angstrom\n`;
  atom_positions.forEach(pos => {
    input += `  ${pos.type}  ${pos.x}  ${pos.y}  ${pos.z}\n`;
  });

  // Add K_POINTS
  input += `\nK_POINTS ${k_points.type}\n`;
  if (k_points.type === 'automatic' && k_points.grid) {
    input += `  ${k_points.grid.join(' ')}\n`;
  }

  return input;
}

// Tool to validate QE parameters
export const validateQEParametersTool = new DynamicStructuredTool({
  name: "validate_qe_parameters",
  description: "Validate Quantum Espresso input parameters and provide recommendations",
  schema: z.object({
    ecutwfc: z.number(),
    ecutrho: z.number().optional(),
    nat: z.number(),
    ntyp: z.number(),
  }),
  func: async ({ ecutwfc, ecutrho, nat, ntyp }) => {
    const warnings = [];
    
    if (ecutwfc < 20) {
      warnings.push("ecutwfc is very low, consider using at least 30 Ry for accurate results");
    }
    
    if (ecutrho && ecutrho < 4 * ecutwfc) {
      warnings.push("ecutrho should typically be at least 4 times ecutwfc for norm-conserving pseudopotentials");
    }
    
    if (nat > 100) {
      warnings.push("Large number of atoms - consider using parallel execution");
    }
    
    return warnings.length > 0 
      ? `Validation warnings: ${warnings.join('; ')}` 
      : "Parameters look good!";
  },
});

// Tool to suggest calculation types
export const suggestCalculationTool = new DynamicStructuredTool({
  name: "suggest_calculation_type",
  description: "Suggest appropriate calculation type based on user requirements",
  schema: z.object({
    purpose: z.string().describe("What the user wants to calculate"),
  }),
  func: async ({ purpose }) => {
    const purposeLower = purpose.toLowerCase();
    
    if (purposeLower.includes('band') || purposeLower.includes('electronic structure')) {
      return "For band structure calculations, use 'scf' first, then 'nscf', then 'bands'";
    } else if (purposeLower.includes('relax') || purposeLower.includes('optimize') || purposeLower.includes('geometry')) {
      if (purposeLower.includes('lattice') || purposeLower.includes('cell') || purposeLower.includes('volume')) {
        return "Use 'vc-relax' for variable-cell relaxation (optimizes both atomic positions and cell parameters)";
      }
      return "Use 'relax' for geometry optimization (fixed cell)";
    } else if (purposeLower.includes('energy') || purposeLower.includes('ground state')) {
      return "Use 'scf' for self-consistent field calculation to get ground state energy";
    }
    
    return "Based on your description, 'scf' is a good starting point. Specify if you need relaxation or band structure.";
  },
});
