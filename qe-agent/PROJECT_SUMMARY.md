# QE Input Generator - Project Summary

## 🎯 Project Overview

An AI-powered web application that uses LangChain agents to generate Quantum ESPRESSO (QE) input files through a conversational interface built with Next.js.

## ✨ Key Features

### 1. **AI Agent with LangChain**
- Uses OpenAI GPT-4 (configurable model)
- Three specialized tools:
  - `generate_qe_input`: Creates properly formatted QE input files
  - `validate_qe_parameters`: Validates parameters and provides recommendations
  - `suggest_calculation_type`: Suggests calculation types based on user needs

### 2. **Natural Language Interface**
- Chat with the AI to describe your calculation needs
- AI asks clarifying questions when needed
- Conversational refinement of parameters

### 3. **Modern Next.js Frontend**
- Beautiful, responsive UI with Tailwind CSS
- Real-time chat interface
- Dark mode support
- File preview and download

### 4. **QE-Specific Features**
- Supports multiple calculation types (scf, relax, vc-relax, nscf, bands)
- Parameter validation with physics-based recommendations
- Proper formatting for PWscf input files
- Examples included for common calculations

## 📁 Project Structure

```
/workspace/qe-agent/
├── app/
│   ├── api/generate/route.ts       # API endpoint for LangChain agent
│   ├── layout.tsx                  # Root layout with metadata
│   ├── page.tsx                    # Main page component
│   └── globals.css                 # Global styles with Tailwind
│
├── components/
│   └── ChatInterface.tsx           # Main chat UI component
│
├── lib/
│   ├── agent.ts                    # LangChain agent configuration
│   └── qe-tools.ts                 # QE-specific tools and schemas
│
├── examples/
│   ├── silicon_scf.in              # Example: Silicon SCF calculation
│   └── water_relax.in              # Example: Water molecule relaxation
│
├── Configuration Files
│   ├── package.json                # Dependencies and scripts
│   ├── tsconfig.json               # TypeScript configuration
│   ├── next.config.js              # Next.js configuration
│   ├── tailwind.config.js          # Tailwind CSS configuration
│   ├── postcss.config.js           # PostCSS configuration
│   └── .eslintrc.json              # ESLint configuration
│
├── Documentation
│   ├── README.md                   # Comprehensive documentation
│   ├── QUICKSTART.md              # Quick start guide
│   └── PROJECT_SUMMARY.md         # This file
│
└── Environment
    ├── .env.example                # Environment variables template
    ├── .env                        # Your API keys (not committed)
    └── .gitignore                  # Git ignore rules
```

## 🛠️ Technology Stack

| Component | Technology | Purpose |
|-----------|-----------|---------|
| Frontend Framework | Next.js 15 | React-based web framework |
| Language | TypeScript | Type-safe JavaScript |
| Styling | Tailwind CSS | Utility-first CSS framework |
| AI Framework | LangChain | Agent orchestration |
| LLM | OpenAI GPT-4o-mini | Natural language understanding |
| Validation | Zod | Schema validation |
| API | Next.js API Routes | Backend API endpoints |

## 🚀 How It Works

### 1. User Interaction Flow

```
User Input → ChatInterface → API Route → LangChain Agent
                                              ↓
                                         QE Tools
                                              ↓
                                    Generated Input File
                                              ↓
                                    Display & Download
```

### 2. LangChain Agent Architecture

The agent uses a ReAct (Reasoning + Acting) pattern:

1. **Reasoning**: Agent analyzes user request
2. **Tool Selection**: Chooses appropriate tool(s)
3. **Action**: Executes tool with parameters
4. **Observation**: Processes tool output
5. **Response**: Generates user-friendly response

### 3. QE Tools

#### generate_qe_input
- **Input**: Complete QE parameters (validated by Zod schema)
- **Output**: Properly formatted PWscf input file
- **Features**: 
  - Handles all namelists (&CONTROL, &SYSTEM, &ELECTRONS, &IONS, &CELL)
  - Formats atomic species and positions
  - Configures k-points sampling

#### validate_qe_parameters
- **Input**: Key parameters (ecutwfc, ecutrho, nat, ntyp)
- **Output**: Validation warnings and recommendations
- **Checks**:
  - Energy cutoff values
  - ecutrho vs ecutwfc ratio
  - System size considerations

#### suggest_calculation_type
- **Input**: User's description of what they want to calculate
- **Output**: Recommended calculation type with explanation
- **Supports**: Energy, relaxation, bands, DOS calculations

## 📝 Example Usage

### Example 1: Simple Request
```
User: "Create an SCF calculation for silicon"

Agent: 
1. Asks about number of atoms, lattice parameters
2. Validates ecutwfc choice
3. Generates complete input file
4. Provides download link
```

### Example 2: Complex Request
```
User: "I want to optimize a water molecule and then calculate its band structure"

Agent:
1. Suggests doing relax first, then scf, then bands
2. Asks about cell size (for molecule in vacuum)
3. Generates relax.in with appropriate parameters
4. Explains next steps for band structure
```

## 🎨 UI Features

- **Chat History**: Maintains conversation context
- **Loading States**: Animated indicators during generation
- **File Preview**: Shows generated input in terminal-style display
- **One-Click Download**: Downloads as `pwscf.in`
- **Clear Button**: Reset conversation
- **Responsive Design**: Works on desktop and mobile
- **Dark Mode**: Automatic theme detection

## 🔧 Configuration Options

### Environment Variables
```bash
OPENAI_API_KEY=sk-...           # Required: Your OpenAI API key
OPENAI_MODEL=gpt-4o-mini        # Optional: Model selection
```

### Customization Points

1. **System Message** (`lib/agent.ts`):
   - Modify agent behavior and expertise
   - Add domain-specific instructions

2. **Tools** (`lib/qe-tools.ts`):
   - Add new QE-specific tools
   - Extend validation logic
   - Add more calculation types

3. **UI** (`components/ChatInterface.tsx`):
   - Customize colors and styling
   - Add new features (history, templates)
   - Modify chat behavior

## 📊 Supported QE Features

### Calculation Types
- ✅ SCF (self-consistent field)
- ✅ Relax (geometry optimization)
- ✅ VC-Relax (variable-cell relaxation)
- ✅ NSCF (non-self-consistent)
- ✅ Bands (band structure)

### Parameters
- ✅ System size (nat, ntyp)
- ✅ Energy cutoffs (ecutwfc, ecutrho)
- ✅ Lattice parameters (ibrav, celldm, a/b/c)
- ✅ Electronic parameters (occupations, smearing)
- ✅ Convergence (mixing_beta, conv_thr)
- ✅ Atomic species (mass, pseudopotentials)
- ✅ Atomic positions (Cartesian, crystal, alat)
- ✅ K-points (automatic, gamma, crystal)

## 🚦 Getting Started

1. **Install dependencies**:
   ```bash
   cd /workspace/qe-agent
   npm install
   ```

2. **Configure API key**:
   ```bash
   echo "OPENAI_API_KEY=sk-your-key" > .env
   ```

3. **Run development server**:
   ```bash
   npm run dev
   ```

4. **Open browser**:
   ```
   http://localhost:3000
   ```

## 🔄 Development Workflow

### Testing Changes
```bash
npm run dev          # Start dev server with hot reload
```

### Building for Production
```bash
npm run build        # Create production build
npm start            # Run production server
```

### Linting
```bash
npm run lint         # Check code quality
```

## 📚 Learning Resources

- **QE Documentation**: https://www.quantum-espresso.org/documentation/
- **LangChain JS**: https://js.langchain.com/docs/
- **Next.js**: https://nextjs.org/docs
- **Tailwind CSS**: https://tailwindcss.com/docs

## 🎯 Future Enhancements

Potential improvements:
- [ ] Support for more QE packages (ph.x, bands.x, etc.)
- [ ] Template library for common materials
- [ ] Parameter optimization suggestions
- [ ] Integration with Materials Project API
- [ ] Batch input file generation
- [ ] Visualization of crystal structures
- [ ] Export to other formats (VASP, ABINIT, etc.)
- [ ] User authentication and saved conversations
- [ ] Pseudopotential recommendations

## 🤝 Contributing

The codebase is well-structured for contributions:
- Add new tools in `lib/qe-tools.ts`
- Extend agent capabilities in `lib/agent.ts`
- Improve UI in `components/ChatInterface.tsx`
- Add examples in `examples/`

## 📄 License

MIT License - Free to use for research and education

## 🙏 Acknowledgments

Built with:
- Next.js team for the amazing framework
- LangChain developers for the agent framework
- OpenAI for GPT models
- Quantum ESPRESSO developers for the DFT package

---

**Project Status**: ✅ Complete and ready to use!

**Last Updated**: October 30, 2025
