import { useState, useEffect } from 'react'

interface ProviderConfig {
  provider: string;
  api_key?: string;
  base_url?: string;
  model: string;
  configured: boolean;
}

export default function Settings() {
  const [configs, setConfigs] = useState<Record<string, ProviderConfig>>({})
  const [activeTab, setActiveTab] = useState('api-keys')

  const providerList = [
    { id: 'openai', name: 'OpenAI', models: ['gpt-4-turbo', 'gpt-3.5-turbo'] },
    { id: 'claude', name: 'Claude', models: ['claude-3-opus', 'claude-3-sonnet'] },
    { id: 'gemini', name: 'Gemini', models: ['gemini-pro', 'gemini-pro-vision'] },
    { id: 'mistral', name: 'Mistral', models: ['mistral-large', 'mistral-medium'] },
    { id: 'deepseek', name: 'DeepSeek', models: ['deepseek-chat'] },
    { id: 'qwen', name: 'Qwen', models: ['qwen-turbo', 'qwen-plus'] },
    { id: 'siliconflow', name: 'SiliconFlow', models: ['Qwen/Qwen2-7B-Instruct'] },
    { id: 'openrouter', name: 'OpenRouter', models: ['anthropic/claude-3-opus', 'openai/gpt-4'] },
    { id: 'ollama', name: 'Ollama', models: ['llama2', 'mistral'] }
  ]

  useEffect(() => {
    providerList.forEach(p => {
      fetch(`/api/ai/config/${p.id}`)
        .then(res => res.json())
        .then(data => setConfigs(prev => ({ ...prev, [p.id]: data })))
        .catch(() => {})
    })
  }, [])

  const handleSave = async (provider: string, config: ProviderConfig) => {
    try {
      await fetch('/api/ai/config', {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(config)
      })
      alert('Configuration saved!')
    } catch (error) {
      console.error('Failed to save config:', error)
    }
  }

  return (
    <div className="page">
      <div className="page-header">
        <h2>Settings</h2>
        <p>Configure AI providers and system preferences</p>
      </div>

      <div className="card">
        <div style={{ display: 'flex', gap: '1rem', marginBottom: '1rem', borderBottom: '1px solid var(--border)', paddingBottom: '1rem' }}>
          <button
            className={`btn ${activeTab === 'api-keys' ? 'btn-primary' : 'btn-secondary'}`}
            onClick={() => setActiveTab('api-keys')}
          >
            API Keys
          </button>
          <button
            className={`btn ${activeTab === 'providers' ? 'btn-primary' : 'btn-secondary'}`}
            onClick={() => setActiveTab('providers')}
          >
            Providers
          </button>
          <button
            className={`btn ${activeTab === 'about' ? 'btn-primary' : 'btn-secondary'}`}
            onClick={() => setActiveTab('about')}
          >
            About
          </button>
        </div>

        {activeTab === 'api-keys' && (
          <div className="provider-grid">
            {providerList.map(p => (
              <div key={p.id} className="provider-card">
                <h4>{p.name}</h4>
                <p className="status">
                  {configs[p.id]?.configured ? 'Configured' : 'Not configured'}
                </p>
                <div className="form-group" style={{ marginTop: '0.5rem' }}>
                  <label>API Key</label>
                  <input
                    type="password"
                    placeholder="sk-..."
                    defaultValue={configs[p.id]?.api_key || ''}
                    onChange={e => setConfigs({
                      ...configs,
                      [p.id]: { ...configs[p.id], provider: p.id, api_key: e.target.value, model: configs[p.id]?.model || p.models[0] }
                    })}
                  />
                </div>
                <button
                  className="btn btn-primary"
                  style={{ width: '100%', marginTop: '0.5rem' }}
                  onClick={() => handleSave(p.id, configs[p.id] || { provider: p.id, model: p.models[0] })}
                >
                  Save
                </button>
              </div>
            ))}
          </div>
        )}

        {activeTab === 'providers' && (
          <div>
            <h3 style={{ marginBottom: '1rem' }}>Available AI Providers</h3>
            <table style={{ width: '100%', borderCollapse: 'collapse' }}>
              <thead>
                <tr style={{ borderBottom: '1px solid var(--border)' }}>
                  <th style={{ textAlign: 'left', padding: '0.5rem' }}>Provider</th>
                  <th style={{ textAlign: 'left', padding: '0.5rem' }}>Status</th>
                  <th style={{ textAlign: 'left', padding: '0.5rem' }}>Models</th>
                </tr>
              </thead>
              <tbody>
                {providerList.map(p => (
                  <tr key={p.id} style={{ borderBottom: '1px solid var(--border)' }}>
                    <td style={{ padding: '0.75rem' }}>{p.name}</td>
                    <td style={{ padding: '0.75rem' }}>
                      <span className={`status-dot ${configs[p.id]?.configured ? '' : 'warning'}`}></span>
                    </td>
                    <td style={{ padding: '0.75rem', color: 'var(--text-secondary)' }}>
                      {p.models.join(', ')}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}

        {activeTab === 'about' && (
          <div>
            <h3>BioDockify v2.3.0</h3>
            <p style={{ marginTop: '1rem', color: 'var(--text-secondary)' }}>
              Production-ready molecular docking software with Discovery Studio-inspired UI.
            </p>
            <ul style={{ marginTop: '1rem', paddingLeft: '1.5rem', color: 'var(--text-secondary)' }}>
              <li>AutoDock Vina for molecular docking</li>
              <li>RDKit for ligand preparation</li>
              <li>OpenMM for molecular dynamics</li>
              <li>Multi-provider AI assistant</li>
            </ul>
          </div>
        )}
      </div>
    </div>
  )
}
