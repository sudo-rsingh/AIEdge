# Quick Start Guide

## 1. Set up your API key

```bash
cd /workspace/qe-agent
echo "OPENAI_API_KEY=your_key_here" > .env
echo "OPENAI_MODEL=gpt-4o-mini" >> .env
```

## 2. Install dependencies (if not already done)

```bash
npm install
```

## 3. Run the development server

```bash
npm run dev
```

## 4. Open your browser

Navigate to [http://localhost:3000](http://localhost:3000)

## 5. Try these example prompts

### Simple SCF calculation
```
Create an SCF calculation for silicon with 2 atoms
```

### Geometry optimization
```
I want to optimize the geometry of a water molecule with 1 oxygen and 2 hydrogen atoms
```

### Band structure
```
Help me set up a band structure calculation for graphene
```

### Custom parameters
```
I need an SCF calculation for:
- Material: Aluminum (FCC)
- 4 atoms in the unit cell
- Energy cutoff: 50 Rydberg
- K-points: 12x12x12 grid
```

## Common Issues

### Missing API Key
If you see errors about missing API key, make sure:
1. `.env` file exists in `/workspace/qe-agent/`
2. File contains `OPENAI_API_KEY=sk-...`
3. Restart the dev server after adding the key

### Port Already in Use
If port 3000 is busy:
```bash
npm run dev -- -p 3001
```

### Dependencies Issues
If you encounter module errors:
```bash
rm -rf node_modules package-lock.json
npm install
```

## Next Steps

1. Explore the example input files in `/workspace/qe-agent/examples/`
2. Customize the agent's system message in `lib/agent.ts`
3. Add more QE-specific tools in `lib/qe-tools.ts`
4. Modify the UI in `components/ChatInterface.tsx`

Enjoy creating QE input files! 🚀
