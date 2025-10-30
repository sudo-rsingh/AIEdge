import { NextRequest, NextResponse } from 'next/server';
import { createQEAgent, QE_SYSTEM_MESSAGE } from '@/lib/agent';
import { HumanMessage, SystemMessage } from '@langchain/core/messages';

export const runtime = 'nodejs';
export const dynamic = 'force-dynamic';

export async function POST(req: NextRequest) {
  try {
    const { message, conversationHistory = [] } = await req.json();

    if (!message) {
      return NextResponse.json(
        { error: 'Message is required' },
        { status: 400 }
      );
    }

    // Check for API key
    if (!process.env.OPENAI_API_KEY) {
      return NextResponse.json(
        { error: 'OpenAI API key not configured' },
        { status: 500 }
      );
    }

    // Create the agent
    const agent = createQEAgent();

    // Build message history
    const messages: any[] = [
      new SystemMessage(QE_SYSTEM_MESSAGE),
    ];

    // Add conversation history
    conversationHistory.forEach((msg: any) => {
      messages.push(new HumanMessage(msg.content));
    });

    // Add current message
    messages.push(new HumanMessage(message));

    // Invoke the agent
    const result = await agent.invoke({
      messages: messages,
    });

    // Extract the final response
    const lastMessage = result.messages[result.messages.length - 1];
    const response = lastMessage.content;

    // Extract generated file if present
    let generatedFile = null;
    const fileMatch = response.match(/```(?:text|pwscf)?\n([\s\S]+?)\n```/);
    if (fileMatch) {
      generatedFile = fileMatch[1];
    }

    return NextResponse.json({
      response: response,
      generatedFile: generatedFile,
      messages: result.messages.map((m: any) => ({
        type: m.constructor.name,
        content: m.content,
      })),
    });

  } catch (error: any) {
    console.error('Error generating QE input:', error);
    return NextResponse.json(
      { 
        error: 'Failed to generate QE input',
        details: error.message 
      },
      { status: 500 }
    );
  }
}
