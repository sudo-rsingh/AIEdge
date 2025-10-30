# QE Input Generator - Architecture Overview

## 🏗️ System Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                         Browser                              │
│  ┌────────────────────────────────────────────────────┐     │
│  │         ChatInterface Component                     │     │
│  │  - User input                                       │     │
│  │  - Message display                                  │     │
│  │  - File preview & download                          │     │
│  └────────────────┬───────────────────────────────────┘     │
└───────────────────┼──────────────────────────────────────────┘
                    │ HTTP POST /api/generate
                    ▼
┌─────────────────────────────────────────────────────────────┐
│                    Next.js API Route                         │
│  ┌────────────────────────────────────────────────────┐     │
│  │         /app/api/generate/route.ts                  │     │
│  │  - Receives user message                            │     │
│  │  - Manages conversation history                     │     │
│  │  - Calls LangChain agent                            │     │
│  │  - Returns response + generated file                │     │
│  └────────────────┬───────────────────────────────────┘     │
└───────────────────┼──────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────┐
│                   LangChain Agent                            │
│  ┌────────────────────────────────────────────────────┐     │
│  │         lib/agent.ts                                │     │
│  │  - ReAct agent (Reasoning + Acting)                │     │
│  │  - OpenAI GPT-4o-mini LLM                          │     │
│  │  - System message with QE expertise                │     │
│  └────────────────┬───────────────────────────────────┘     │
└───────────────────┼──────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────┐
│                      QE Tools                                │
│  ┌────────────────────────────────────────────────────┐     │
│  │         lib/qe-tools.ts                             │     │
│  │                                                      │     │
│  │  Tool 1: generate_qe_input                          │     │
│  │  - Input: QE parameters (Zod validated)             │     │
│  │  - Output: Formatted PWscf input file               │     │
│  │                                                      │     │
│  │  Tool 2: validate_qe_parameters                     │     │
│  │  - Input: ecutwfc, ecutrho, nat, ntyp               │     │
│  │  - Output: Validation warnings                      │     │
│  │                                                      │     │
│  │  Tool 3: suggest_calculation_type                   │     │
│  │  - Input: User's description                        │     │
│  │  - Output: Recommended calc type                    │     │
│  └─────────────────────────────────────────────────────┘     │
└─────────────────────────────────────────────────────────────┘
```

## 📊 Data Flow

### 1. User Request Flow
```
User: "Create SCF calculation for silicon"
  ↓
ChatInterface.tsx
  ↓ (POST /api/generate)
API Route (route.ts)
  ↓ (invoke agent with message)
LangChain Agent
  ↓ (reasoning: need to generate QE input)
  ↓ (tool selection: generate_qe_input)
QE Tools (qe-tools.ts)
  ↓ (execute tool)
  ↓ (return formatted input file)
Agent
  ↓ (format response for user)
API Route
  ↓ (JSON response)
ChatInterface
  ↓ (display message + file)
User sees: Response + Download button
```

### 2. Agent Reasoning Pattern (ReAct)

```
┌─────────────────────────────────────────────┐
│ 1. THOUGHT                                  │
│    "User wants SCF for silicon, I need      │
│     atomic positions and lattice params"    │
└────────────────┬────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────┐
│ 2. ACTION                                   │
│    suggest_calculation_type("scf")          │
│    → Returns: "SCF is correct choice"       │
└────────────────┬────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────┐
│ 3. THOUGHT                                  │
│    "Need ecutwfc, I'll suggest 30 Ry"       │
└────────────────┬────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────┐
│ 4. ACTION                                   │
│    validate_qe_parameters(ecutwfc=30, ...)  │
│    → Returns: "Parameters look good"        │
└────────────────┬────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────┐
│ 5. ACTION                                   │
│    generate_qe_input({calculation: 'scf',   │
│                       nat: 2, ...})         │
│    → Returns: Complete input file           │
└────────────────┬────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────┐
│ 6. ANSWER                                   │
│    "Here's your SCF input file for silicon" │
└─────────────────────────────────────────────┘
```

## 🔧 Component Breakdown

### Frontend Layer

**ChatInterface.tsx**
- **State Management**:
  - `messages`: Conversation history
  - `input`: Current user input
  - `loading`: Request in progress
  - `generatedFile`: Latest generated QE input

- **Functions**:
  - `handleSubmit()`: Send message to API
  - `handleDownload()`: Download generated file
  - `handleClear()`: Reset conversation
  - `scrollToBottom()`: Auto-scroll to latest

- **UI Components**:
  - Message bubbles (user/assistant)
  - Input form
  - Loading animation
  - File preview panel
  - Download button

### Backend Layer

**API Route (app/api/generate/route.ts)**
- **Responsibilities**:
  - Validate request
  - Check API key
  - Build message history
  - Invoke agent
  - Extract generated file
  - Return response

- **Request Format**:
  ```json
  {
    "message": "Create SCF for silicon",
    "conversationHistory": [...]
  }
  ```

- **Response Format**:
  ```json
  {
    "response": "Here's your input file...",
    "generatedFile": "&CONTROL\n...",
    "messages": [...]
  }
  ```

### Agent Layer

**Agent (lib/agent.ts)**
- **Configuration**:
  - Model: gpt-4o-mini (configurable)
  - Temperature: 0.1 (low for consistency)
  - Tools: 3 QE-specific tools

- **System Message**:
  - Defines agent role (QE expert)
  - Lists available tools
  - Provides calculation type info
  - Sets helpful behavior

### Tools Layer

**QE Tools (lib/qe-tools.ts)**

1. **generate_qe_input**
   - **Schema**: 20+ validated fields
   - **Logic**:
     ```
     Build &CONTROL namelist
     Build &SYSTEM namelist
     Build &ELECTRONS namelist
     Add &IONS if relax
     Add &CELL if vc-relax
     Format ATOMIC_SPECIES
     Format ATOMIC_POSITIONS
     Format K_POINTS
     Return complete file
     ```

2. **validate_qe_parameters**
   - **Checks**:
     - ecutwfc >= 20 Ry
     - ecutrho >= 4 * ecutwfc
     - Large system warnings
   - **Returns**: Warning messages or "OK"

3. **suggest_calculation_type**
   - **Pattern Matching**:
     - "band" → scf + nscf + bands
     - "relax" → relax or vc-relax
     - "energy" → scf
   - **Returns**: Recommendation + explanation

## 🗂️ File Organization

```
qe-agent/
│
├── Frontend (User Interface)
│   ├── app/
│   │   ├── page.tsx              # Main page
│   │   ├── layout.tsx            # App layout
│   │   └── globals.css           # Global styles
│   └── components/
│       └── ChatInterface.tsx     # Chat UI
│
├── Backend (API & Logic)
│   ├── app/api/generate/
│   │   └── route.ts              # API endpoint
│   └── lib/
│       ├── agent.ts              # LangChain agent
│       └── qe-tools.ts           # QE tools
│
├── Configuration
│   ├── package.json              # Dependencies
│   ├── tsconfig.json             # TypeScript
│   ├── next.config.js            # Next.js
│   ├── tailwind.config.js        # Tailwind
│   └── .env                      # Environment
│
├── Documentation
│   ├── README.md                 # Full docs
│   ├── QUICKSTART.md             # Quick start
│   ├── PROJECT_SUMMARY.md        # Overview
│   └── ARCHITECTURE.md           # This file
│
└── Examples
    ├── silicon_scf.in            # SCF example
    └── water_relax.in            # Relax example
```

## 🔐 Security & Best Practices

### Environment Variables
- ✅ API keys in `.env` (not committed)
- ✅ `.env.example` provided as template
- ✅ Server-side only (Next.js API routes)

### Input Validation
- ✅ Zod schemas for type safety
- ✅ Parameter range checking
- ✅ Physics-based validation

### Error Handling
- ✅ Try-catch blocks
- ✅ User-friendly error messages
- ✅ Graceful degradation

### Code Quality
- ✅ TypeScript for type safety
- ✅ ESLint for code quality
- ✅ Consistent formatting

## 🎨 Styling Architecture

### Tailwind CSS Classes
- **Colors**: `bg-blue-600`, `text-white`, etc.
- **Layout**: `flex`, `grid`, `max-w-6xl`, etc.
- **Responsive**: Mobile-first design
- **Dark Mode**: `dark:` variants

### Custom Styles
- **Scrollbar**: Custom thin scrollbar
- **Animations**: Loading dots with stagger
- **Gradients**: Header and background

## 🔄 State Management

### Client State (React)
```typescript
const [messages, setMessages] = useState<Message[]>([])
const [input, setInput] = useState('')
const [loading, setLoading] = useState(false)
const [generatedFile, setGeneratedFile] = useState<string | null>(null)
```

### Server State (API)
- Conversation history passed in request
- Stateless API (no session storage)
- Each request is independent

## 📈 Performance Optimizations

- **Streaming**: Could add streaming responses
- **Caching**: Could cache common queries
- **Lazy Loading**: Components load on demand
- **Code Splitting**: Next.js automatic splitting

## 🧪 Testing Strategy (Future)

### Unit Tests
- QE tool functions
- Parameter validation
- File generation logic

### Integration Tests
- API endpoints
- Agent tool calling
- End-to-end flows

### E2E Tests
- User interaction flows
- File download
- Error handling

## 🚀 Deployment Options

### Vercel (Recommended)
```bash
vercel deploy
```

### Docker
```dockerfile
FROM node:18-alpine
WORKDIR /app
COPY . .
RUN npm install
RUN npm run build
CMD ["npm", "start"]
```

### Traditional Server
```bash
npm run build
npm start
```

---

**Architecture Status**: Production-ready ✅

This architecture is:
- ✅ Scalable (stateless API)
- ✅ Maintainable (modular design)
- ✅ Extensible (easy to add tools)
- ✅ Type-safe (TypeScript + Zod)
- ✅ User-friendly (modern UI)
