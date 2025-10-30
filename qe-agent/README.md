# QE Input Generator - AI Agent

An AI-powered web application for generating Quantum ESPRESSO (QE) input files using LangChain agents and Next.js.

## Features

- 🤖 **AI-Powered Generation**: Uses LangChain agents with GPT-4 to understand natural language requests
- 💬 **Conversational Interface**: Chat with the AI to refine your QE input parameters
- 📝 **Smart Validation**: Automatic parameter validation and recommendations
- 📥 **Easy Download**: Download generated input files with one click
- 🎨 **Modern UI**: Beautiful, responsive interface built with Next.js and Tailwind CSS
- 🔧 **Multiple Calculation Types**: Supports SCF, relax, vc-relax, nscf, and bands calculations

## Prerequisites

- Node.js 18+ 
- npm or yarn
- OpenAI API key

## Installation

1. Clone the repository or navigate to the project directory:
```bash
cd qe-agent
```

2. Install dependencies:
```bash
npm install
```

3. Set up environment variables:
```bash
cp .env.example .env
```

4. Edit `.env` and add your OpenAI API key:
```
OPENAI_API_KEY=your_openai_api_key_here
OPENAI_MODEL=gpt-4o-mini
```

## Usage

### Development Mode

Run the development server:
```bash
npm run dev
```

Open [http://localhost:3000](http://localhost:3000) in your browser.

### Production Build

Build and run the production version:
```bash
npm run build
npm start
```

## How to Use

1. **Start a conversation**: The AI assistant will greet you and ask what calculation you want to perform.

2. **Describe your calculation**: Use natural language to describe your needs. Examples:
   - "I want to run an SCF calculation for silicon with 2 atoms"
   - "Create a geometry optimization for a water molecule"
   - "I need to calculate the band structure of graphene"

3. **Provide details**: The AI will ask for necessary parameters like:
   - Number of atoms and atom types
   - Lattice parameters
   - Energy cutoffs
   - K-point sampling
   - Pseudopotentials

4. **Review and download**: Once generated, review the input file and download it using the "Download" button.

## Example Conversations

### Example 1: Simple SCF Calculation
```
User: Create an SCF calculation for silicon with 2 atoms

AI: I'll help you create an SCF calculation for silicon. Let me gather the necessary information...
[The AI will ask about lattice parameters, pseudopotentials, and other details]
```

### Example 2: Geometry Optimization
```
User: I need to optimize the geometry of a water molecule

AI: I'll set up a relaxation calculation for H2O. For molecules, we typically use:
- ibrav = 0 (free lattice)
- A large cubic cell to avoid periodic interactions
[Continues to gather information and generate the input]
```

## Project Structure

```
qe-agent/
├── app/
│   ├── api/
│   │   └── generate/
│   │       └── route.ts          # API endpoint for agent
│   ├── layout.tsx                # Root layout
│   ├── page.tsx                  # Home page
│   └── globals.css               # Global styles
├── components/
│   └── ChatInterface.tsx         # Main chat component
├── lib/
│   ├── agent.ts                  # LangChain agent setup
│   └── qe-tools.ts              # QE-specific tools and schemas
├── .env.example                  # Environment variables template
├── next.config.js               # Next.js configuration
├── tailwind.config.js           # Tailwind CSS configuration
└── tsconfig.json                # TypeScript configuration
```

## LangChain Tools

The agent has access to three specialized tools:

1. **generate_qe_input**: Creates properly formatted QE input files
2. **validate_qe_parameters**: Validates parameters and provides recommendations
3. **suggest_calculation_type**: Suggests appropriate calculation types based on user needs

## Supported Calculation Types

- **scf**: Self-consistent field calculation (ground state energy)
- **relax**: Geometry optimization with fixed cell parameters
- **vc-relax**: Variable-cell relaxation (optimizes both atoms and cell)
- **nscf**: Non-self-consistent calculation (for DOS, bands)
- **bands**: Band structure calculation

## Common QE Parameters

### Required Parameters
- `calculation`: Type of calculation (scf, relax, etc.)
- `nat`: Number of atoms
- `ntyp`: Number of atomic species
- `ecutwfc`: Kinetic energy cutoff for wavefunctions (Ry)
- `ibrav`: Bravais lattice index

### Typical Values
- `ecutwfc`: 30-80 Ry (depends on pseudopotentials)
- `ecutrho`: 4-12 times ecutwfc
- `mixing_beta`: 0.3-0.7 (for convergence)
- `conv_thr`: 1.0e-6 to 1.0e-8 (convergence threshold)

## Troubleshooting

### API Key Error
If you see "OpenAI API key not configured":
- Check that `.env` file exists in the project root
- Verify `OPENAI_API_KEY` is set correctly
- Restart the development server

### Module Not Found
If you get module errors:
```bash
rm -rf node_modules package-lock.json
npm install
```

### Build Errors
For TypeScript errors during build:
```bash
npm run build
```
Check the error messages and ensure all types are properly defined.

## Technologies Used

- **Next.js 15**: React framework for production
- **TypeScript**: Type-safe JavaScript
- **Tailwind CSS**: Utility-first CSS framework
- **LangChain**: Framework for building LLM applications
- **OpenAI GPT-4**: Language model for understanding and generation
- **Zod**: TypeScript-first schema validation

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

MIT License - feel free to use this project for your research or educational purposes.

## Acknowledgments

- Quantum ESPRESSO developers for the excellent DFT package
- LangChain team for the agent framework
- OpenAI for GPT models

## Further Reading

- [Quantum ESPRESSO Documentation](https://www.quantum-espresso.org/documentation/)
- [LangChain Documentation](https://js.langchain.com/docs/)
- [Next.js Documentation](https://nextjs.org/docs)

## Support

For issues related to:
- **QE Input Generation**: Open an issue in this repository
- **Quantum ESPRESSO**: Visit [QE forums](https://www.quantum-espresso.org/forum/)
- **LangChain**: Check [LangChain GitHub](https://github.com/langchain-ai/langchainjs)
