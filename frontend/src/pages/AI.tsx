import { useState, useEffect } from 'react'

interface Message {
  role: 'user' | 'ai';
  content: string;
}

interface Provider {
  id: string;
  name: string;
}

export default function AIAssistant() {
  const [messages, setMessages] = useState<Message[]>([])
  const [input, setInput] = useState('')
  const [selectedProvider, setSelectedProvider] = useState('openai')
  const [providers, setProviders] = useState<Provider[]>([])
  const [isLoading, setIsLoading] = useState(false)

  useEffect(() => {
    fetch('/api/ai/providers')
      .then(res => res.json())
      .then(data => setProviders(data.providers || []))
      .catch(console.error)
  }, [])

  const handleSend = async () => {
    if (!input.trim()) return

    const userMessage: Message = { role: 'user', content: input }
    setMessages([...messages, userMessage])
    setInput('')
    setIsLoading(true)

    try {
      const response = await fetch('/api/ai/chat', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          message: input,
          provider: selectedProvider
        })
      })

      const data = await response.json()
      const aiMessage: Message = { role: 'ai', content: data.response }
      setMessages(prev => [...prev, aiMessage])
    } catch (error) {
      console.error('Failed to get AI response:', error)
    }

    setIsLoading(false)
  }

  return (
    <div className="page">
      <div className="page-header">
        <h2>BioDockify AI</h2>
        <p>Ask questions about molecular docking, simulations, and more</p>
      </div>

      <div className="card">
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
          <h3>Chat</h3>
          <select
            value={selectedProvider}
            onChange={e => setSelectedProvider(e.target.value)}
            style={{ padding: '0.5rem', background: 'var(--bg-primary)', color: 'var(--text-primary)', border: '1px solid var(--border)', borderRadius: '4px' }}
          >
            {providers.map(p => (
              <option key={p.id} value={p.id}>{p.name}</option>
            ))}
          </select>
        </div>

        <div className="ai-chat">
          <div className="chat-messages">
            {messages.length === 0 ? (
              <p style={{ color: 'var(--text-secondary)', textAlign: 'center', marginTop: '2rem' }}>
                Ask me anything about molecular docking, MD simulations, or drug discovery!
              </p>
            ) : (
              messages.map((msg, i) => (
                <div key={i} className={`chat-message ${msg.role}`}>
                  <strong>{msg.role === 'user' ? 'You' : 'BioDockify AI'}:</strong> {msg.content}
                </div>
              ))
            )}
            {isLoading && (
              <div className="chat-message ai">
                <em>Thinking...</em>
              </div>
            )}
          </div>

          <div className="chat-input">
            <input
              type="text"
              placeholder="Ask BioDockify AI..."
              value={input}
              onChange={e => setInput(e.target.value)}
              onKeyPress={e => e.key === 'Enter' && handleSend()}
            />
            <button className="btn btn-primary" onClick={handleSend} disabled={isLoading}>
              Send
            </button>
          </div>
        </div>
      </div>
    </div>
  )
}
