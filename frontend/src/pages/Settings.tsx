import { useState, useEffect } from 'react'
import { useTheme } from '@/contexts/ThemeContext'
import { useAccessibility } from '@/contexts/AccessibilityContext'
import { locales, getLocale, setLocale, Locale } from '@/i18n'

const APP_VERSION = '2.4.0'

// Nanobot Plugins (channels/providers)
const PLUGIN_CHANNELS = [
  { id: 'pubchem', name: 'PubChem', icon: '🧪', enabled: true, description: 'Fetch compounds and properties from PubChem' },
  { id: 'pdb', name: 'PDB', icon: '🧬', enabled: true, description: 'Download protein structures from RCSB PDB' },
  { id: 'chembl', name: 'ChEMBL', icon: '💊', enabled: false, description: 'Search ChEMBL database for bioactivity data' },
  { id: 'zinc', name: 'ZINC', icon: '🛒', enabled: false, description: 'Access ZINC compound library for virtual screening' },
  { id: 'pax', name: 'PAX Healthcare', icon: '🏥', enabled: false, description: 'Clinical trial and drug interaction data' },
]

const AI_PROVIDERS = [
  { value: 'ollama', label: 'Ollama (Local)', icon: '🏠', baseUrl: 'http://localhost:11434/v1', needsApiKey: false },
  { value: 'openai', label: 'OpenAI', icon: '🤖', baseUrl: 'https://api.openai.com/v1', needsApiKey: true },
  { value: 'anthropic', label: 'Anthropic Claude', icon: '🧠', baseUrl: 'https://api.anthropic.com/v1', needsApiKey: true },
  { value: 'gemini', label: 'Google Gemini', icon: '✨', baseUrl: 'https://generativelanguage.googleapis.com/v1beta', needsApiKey: true },
  { value: 'deepseek', label: 'DeepSeek', icon: '🔮', baseUrl: 'https://api.deepseek.com/v1', needsApiKey: true },
  { value: 'mistral', label: 'Mistral AI', icon: '🌬️', baseUrl: 'https://api.mistral.ai/v1', needsApiKey: true },
  { value: 'groq', label: 'Groq', icon: '⚡', baseUrl: 'https://api.groq.com/openai/v1', needsApiKey: true },
  { value: 'openrouter', label: 'OpenRouter', icon: '🛣️', baseUrl: 'https://openrouter.ai/api/v1', needsApiKey: true },
  { value: 'siliconflow', label: 'SiliconFlow', icon: '💎', baseUrl: 'https://api.siliconflow.cn/v1', needsApiKey: true },
  { value: 'qwen', label: 'Qwen (Alibaba)', icon: '🌏', baseUrl: 'https://dashscope.aliyuncs.com/compatible-mode/v1', needsApiKey: true },
  { value: 'custom', label: 'Custom (OpenAI-compatible)', icon: '⚙️', baseUrl: '', needsApiKey: true },
]

export function Settings() {
  const { theme, setTheme } = useTheme()
  const isDark = theme === 'dark'
  
  const [activeTab, setActiveTab] = useState<'llm' | 'notifications' | 'plugins' | 'system' | 'about' | 'language' | 'accessibility'>('llm')
  const [showApiKey, setShowApiKey] = useState(false)
  const [saving, setSaving] = useState(false)
  const [testing, setTesting] = useState(false)
  const [autoDetecting, setAutoDetecting] = useState(false)
  const [fetchingModels, setFetchingModels] = useState(false)
  const [testResult, setTestResult] = useState<{ status: string; response?: string; error?: string } | null>(null)
  const [message, setMessage] = useState('')
  const { highContrast, reducedMotion, fontSize, toggleHighContrast, toggleReducedMotion, setFontSize } = useAccessibility()
  const [currentLocale, setCurrentLocale] = useState<Locale>(getLocale())
  const [availableOllamaModels, setAvailableOllamaModels] = useState<string[]>([])

  const [emailConfig, setEmailConfig] = useState({
    enabled: false,
    smtpHost: '',
    smtpPort: '587',
    from: '',
    password: '',
    to: '',
  })

  const [plugins, setPlugins] = useState(PLUGIN_CHANNELS)

  const [llmConfig, setLlmConfig] = useState({
    provider: 'ollama',
    model: '',
    apiKey: '',
    baseUrl: '',
    temperature: '0.7',
    maxTokens: '4096',
  })

  useEffect(() => {
    fetchLLMSettings()
    if (llmConfig.provider === 'ollama') {
      fetchOllamaModels()
    }
  }, [])

  const fetchLLMSettings = async () => {
    try {
      const res = await fetch('/llm/settings')
      const data = await res.json()
      setLlmConfig({
        provider: data.provider || 'ollama',
        model: data.model || '',
        apiKey: data.api_key || '',
        baseUrl: data.base_url || 'http://localhost:11434/v1',
        temperature: String(data.temperature || 0.7),
        maxTokens: String(data.max_tokens || 4096),
      })
    } catch (e) {
      console.error('Failed to fetch LLM settings:', e)
    }
  }

  const fetchOllamaModels = async () => {
    setFetchingModels(true)
    try {
      const res = await fetch('/llm/ollama/models')
      if (res.ok) {
        const data = await res.json()
        const models = data.models || []
        setAvailableOllamaModels(models)
        if (models.length > 0 && !models.includes(llmConfig.model)) {
          setLlmConfig(prev => ({ ...prev, model: models[0] }))
        }
      }
    } catch (e) {
      console.error('Failed to fetch Ollama models:', e)
    } finally {
      setFetchingModels(false)
    }
  }

  const handleProviderChange = (provider: string) => {
    const providerInfo = AI_PROVIDERS.find(p => p.value === provider)
    
    setLlmConfig(prev => ({
      ...prev,
      provider,
      baseUrl: providerInfo?.baseUrl || '',
      apiKey: providerInfo?.needsApiKey ? prev.apiKey : '',
    }))
    setTestResult(null)
    
    if (provider === 'ollama') {
      fetchOllamaModels()
    }
  }

  const handleSave = async () => {
    setSaving(true)
    setMessage('')
    try {
      const res = await fetch('/llm/settings', {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          provider: llmConfig.provider,
          model: llmConfig.model,
          api_key: llmConfig.apiKey,
          base_url: llmConfig.baseUrl,
          temperature: parseFloat(llmConfig.temperature) || 0.7,
          max_tokens: parseInt(llmConfig.maxTokens) || 4096,
        }),
      })
      if (res.ok) {
        setMessage('Settings saved successfully')
      } else {
        setMessage('Failed to save settings')
      }
    } catch {
      setMessage('Failed to save settings')
    } finally {
      setSaving(false)
    }
  }

  const handleTest = async () => {
    setTesting(true)
    setTestResult(null)
    setMessage('')
    try {
      const res = await fetch('/llm/test', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          provider: llmConfig.provider,
          model: llmConfig.model,
          api_key: llmConfig.apiKey,
          base_url: llmConfig.baseUrl,
        }),
      })
      const data = await res.json()
      setTestResult(data)
      if (data.status === 'ok') {
        setMessage('Connection successful!')
      } else {
        setMessage(data.error || 'Connection failed')
      }
    } catch (err: any) {
      setTestResult({ status: 'error', error: err?.message || 'Connection failed' })
      setMessage('Connection failed')
    } finally {
      setTesting(false)
    }
  }

  const handleAutoDetect = async () => {
    setAutoDetecting(true)
    setMessage('')
    try {
      const res = await fetch('/llm/auto-detect')
      const data = await res.json()
      if (data.provider) {
        setLlmConfig({
          provider: data.provider,
          model: data.model,
          apiKey: '',
          baseUrl: data.base_url || '',
          temperature: '0.7',
          maxTokens: '4096',
        })
        setMessage(`Auto-detected: ${data.provider} (${data.model})`)
        if (data.provider === 'ollama') {
          fetchOllamaModels()
        }
      } else {
        setMessage('No LLM provider detected. Install Ollama or add an API key.')
      }
    } catch (e: any) {
      setMessage('Auto-detect failed: ' + e.message)
    } finally {
      setAutoDetecting(false)
    }
  }

  const inputClass = (hasError = false) => `
    w-full px-3 py-2.5 rounded-lg border text-sm transition-all
    ${isDark 
      ? 'bg-gray-700 border-gray-600 text-white placeholder-gray-400 focus:ring-blue-500' 
      : 'bg-white border-gray-300 text-gray-900 placeholder-gray-400 focus:ring-blue-500'
    }
    focus:outline-none focus:ring-2 focus:border-transparent
    ${hasError ? 'border-red-500' : ''}
  `

  const labelClass = `block text-sm font-medium mb-1.5 ${isDark ? 'text-gray-300' : 'text-gray-700'}`
  const hintClass = `text-xs mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`

  const providerInfo = AI_PROVIDERS.find(p => p.value === llmConfig.provider)

  return (
    <div className={`min-h-screen ${isDark ? 'bg-gray-900 text-white' : 'bg-gray-50 text-gray-900'}`}>
      <div className="max-w-5xl mx-auto p-6">
        <div className="flex items-center justify-between mb-8">
          <div>
            <h1 className="text-2xl font-bold">Settings</h1>
            <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
              Configure AI provider and application preferences
            </p>
          </div>
          <span className={`px-3 py-1 rounded-full text-sm border ${
            isDark ? 'bg-gray-800 border-gray-700 text-gray-400' : 'bg-white border-gray-200 text-gray-500'
          }`}>
            v{APP_VERSION}
          </span>
        </div>

        {message && (
          <div className={`mb-4 p-3 rounded-lg text-sm ${
            message.includes('successful') || message.includes('successful')
              ? isDark ? 'bg-green-900/30 text-green-400 border border-green-800' : 'bg-green-50 text-green-700 border border-green-200'
              : isDark ? 'bg-red-900/30 text-red-400 border border-red-800' : 'bg-red-50 text-red-700 border border-red-200'
          }`}>
            {message}
          </div>
        )}

        <div className="flex gap-1 mb-6 p-1 rounded-lg w-fit ${isDark ? 'bg-gray-800' : 'bg-gray-200'} overflow-x-auto">
          {([
            { key: 'llm', label: 'AI Provider', icon: '🤖' },
            { key: 'notifications', label: 'Email', icon: '📧' },
            { key: 'plugins', label: 'Plugins', icon: '🔌' },
            { key: 'system', label: 'System', icon: '⚙️' },
            { key: 'about', label: 'About', icon: 'ℹ️' },
            { key: 'language', label: 'Language', icon: '🌐' },
            { key: 'accessibility', label: 'Accessibility', icon: '♿' },
          ] as const).map(tab => (
            <button
              key={tab.key}
              onClick={() => setActiveTab(tab.key)}
              className={`px-4 py-2 rounded-md text-sm font-medium transition-all whitespace-nowrap ${
                activeTab === tab.key
                  ? isDark ? 'bg-blue-600 text-white' : 'bg-white text-blue-600 shadow-sm'
                  : isDark ? 'text-gray-400 hover:text-white' : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              <span className="mr-1">{tab.icon}</span>
              {tab.label}
            </button>
          ))}
        </div>

        {activeTab === 'llm' && (
          <div className="space-y-6">
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">AI Provider Configuration</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                  Select your preferred LLM provider and model
                </p>
              </div>
              
              <div className="p-6 space-y-5">
                <div>
                  <label className={labelClass}>Provider</label>
                  <div className="grid grid-cols-2 sm:grid-cols-3 md:grid-cols-4 gap-2">
                    {AI_PROVIDERS.map(p => (
                      <button
                        key={p.value}
                        onClick={() => handleProviderChange(p.value)}
                        className={`p-3 rounded-lg border text-left transition-all ${
                          llmConfig.provider === p.value
                            ? isDark 
                              ? 'bg-blue-600/20 border-blue-500 text-blue-400' 
                              : 'bg-blue-50 border-blue-500 text-blue-700'
                            : isDark
                              ? 'bg-gray-700 border-gray-600 text-gray-300 hover:border-gray-500'
                              : 'bg-gray-50 border-gray-200 text-gray-700 hover:border-gray-300'
                        }`}
                      >
                        <div className="flex items-center gap-2">
                          <span>{p.icon}</span>
                          <span className="text-sm font-medium">{p.label}</span>
                        </div>
                      </button>
                    ))}
                  </div>
                </div>

                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <div>
                    <label className={labelClass}>Model</label>
                    <input
                      type="text"
                      value={llmConfig.model}
                      onChange={(e) => setLlmConfig({ ...llmConfig, model: e.target.value })}
                      className={inputClass()}
                      placeholder="e.g. llama3.2, qwen2.5, gpt-4o"
                    />
                    {llmConfig.provider === 'ollama' && (
                      <div className="flex items-center gap-2 mt-1">
                        <button
                          onClick={fetchOllamaModels}
                          disabled={fetchingModels}
                          className={`text-xs px-2 py-1 rounded ${isDark ? 'bg-gray-700 hover:bg-gray-600' : 'bg-gray-100 hover:bg-gray-200'}`}
                        >
                          {fetchingModels ? 'Checking...' : 'Check installed models'}
                        </button>
                        {availableOllamaModels.length > 0 && (
                          <span className="text-xs text-gray-400">
                            Installed: {availableOllamaModels.join(', ')}
                          </span>
                        )}
                      </div>
                    )}
                  </div>

                  <div>
                    <label className={labelClass}>Max Tokens</label>
                    <input
                      type="number"
                      value={llmConfig.maxTokens}
                      onChange={(e) => setLlmConfig({ ...llmConfig, maxTokens: e.target.value })}
                      min={100}
                      max={128000}
                      className={inputClass()}
                      placeholder="4096"
                    />
                    <p className={hintClass}>Maximum response length</p>
                  </div>
                </div>

                <div>
                  <label className={labelClass}>Base URL</label>
                  <input
                    type="text"
                    value={llmConfig.baseUrl}
                    onChange={(e) => setLlmConfig({ ...llmConfig, baseUrl: e.target.value })}
                    className={inputClass()}
                    placeholder={providerInfo?.baseUrl || 'https://api.example.com/v1'}
                  />
                  <p className={hintClass}>
                    {providerInfo?.needsApiKey 
                      ? 'API endpoint for your provider (OpenAI-compatible)'
                      : 'Ollama server URL (default: http://localhost:11434/v1)'}
                  </p>
                </div>

                <div>
                  <label className={labelClass}>API Key</label>
                  <div className="relative">
                    <input
                      type={showApiKey ? 'text' : 'password'}
                      value={llmConfig.apiKey}
                      onChange={(e) => setLlmConfig({ ...llmConfig, apiKey: e.target.value })}
                      className={inputClass() + ' pr-20'}
                      placeholder={providerInfo?.needsApiKey ? 'sk-...' : 'Not required for Ollama'}
                      disabled={!providerInfo?.needsApiKey}
                    />
                    <button
                      type="button"
                      onClick={() => setShowApiKey(!showApiKey)}
                      className={`absolute right-3 top-1/2 -translate-y-1/2 text-xs font-medium ${
                        isDark ? 'text-gray-400 hover:text-white' : 'text-gray-500 hover:text-gray-700'
                      }`}
                    >
                      {showApiKey ? 'Hide' : 'Show'}
                    </button>
                  </div>
                  <p className={hintClass}>
                    {providerInfo?.needsApiKey 
                      ? 'Your API key is stored securely in memory only'
                      : 'No API key required for local Ollama models'}
                  </p>
                </div>

                <div>
                  <label className={labelClass}>Temperature: {llmConfig.temperature}</label>
                  <input
                    type="range"
                    min={0}
                    max={2}
                    step={0.1}
                    value={llmConfig.temperature}
                    onChange={(e) => setLlmConfig({ ...llmConfig, temperature: e.target.value })}
                    className="w-full h-2 rounded-lg appearance-none cursor-pointer bg-gray-300 dark:bg-gray-600"
                  />
                  <div className="flex justify-between text-xs text-gray-400 mt-1">
                    <span>Precise</span>
                    <span>Balanced</span>
                    <span>Creative</span>
                  </div>
                </div>

                <div className="flex gap-3 pt-2">
                  <button
                    onClick={handleAutoDetect}
                    disabled={autoDetecting}
                    className={`px-5 py-2.5 font-medium rounded-lg transition-colors ${
                      isDark 
                        ? 'bg-purple-600 hover:bg-purple-700 disabled:bg-purple-400 text-white' 
                        : 'bg-purple-600 hover:bg-purple-700 disabled:bg-purple-400 text-white'
                    }`}
                  >
                    {autoDetecting ? '⏳ Detecting...' : '🔍 Auto-Detect'}
                  </button>
                  <button
                    onClick={handleSave}
                    disabled={saving}
                    className="px-5 py-2.5 bg-blue-600 hover:bg-blue-700 disabled:bg-blue-400 text-white font-medium rounded-lg transition-colors"
                  >
                    {saving ? 'Saving...' : 'Save Settings'}
                  </button>
                  <button
                    onClick={handleTest}
                    disabled={testing}
                    className={`px-5 py-2.5 font-medium rounded-lg transition-colors ${
                      isDark 
                        ? 'bg-gray-700 hover:bg-gray-600 text-white border border-gray-600' 
                        : 'bg-white hover:bg-gray-50 text-gray-700 border border-gray-300'
                    }`}
                  >
                    {testing ? 'Testing...' : 'Test Connection'}
                  </button>
                </div>

                {testResult && (
                  <div className={`p-4 rounded-lg border ${
                    testResult.status === 'ok'
                      ? isDark ? 'bg-green-900/20 border-green-800 text-green-400' : 'bg-green-50 border-green-200 text-green-700'
                      : isDark ? 'bg-red-900/20 border-red-800 text-red-400' : 'bg-red-50 border-red-200 text-red-700'
                  }`}>
                    <p className="font-semibold text-sm flex items-center gap-2">
                      {testResult.status === 'ok' ? (
                        <>
                          <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                          </svg>
                          Connection successful!
                        </>
                      ) : (
                        <>
                          <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                          </svg>
                          Connection failed
                        </>
                      )}
                    </p>
                    {testResult.response && (
                      <p className="text-sm mt-2 opacity-80">Model: "{testResult.response}"</p>
                    )}
                    {testResult.error && (
                      <p className="text-xs mt-2 font-mono opacity-80">Error: {testResult.error}</p>
                    )}
                  </div>
                )}
              </div>
            </div>
          </div>
        )}

        {activeTab === 'system' && (
          <div className="space-y-6">
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">Theme</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                  Choose your preferred appearance
                </p>
              </div>
              <div className="p-6">
                <div className="flex gap-3">
                  {(['light', 'dark'] as const).map(t => (
                    <button
                      key={t}
                      onClick={() => setTheme(t)}
                      className={`flex-1 p-4 rounded-xl border-2 transition-all ${
                        theme === t
                          ? isDark
                            ? 'border-blue-500 bg-blue-900/30'
                            : 'border-blue-500 bg-blue-50'
                          : isDark
                            ? 'border-gray-700 bg-gray-700/50 hover:border-gray-600'
                            : 'border-gray-200 bg-gray-50 hover:border-gray-300'
                      }`}
                    >
                      <div className="text-3xl mb-2">{t === 'light' ? '☀️' : '🌙'}</div>
                      <div className="font-medium capitalize">{t === 'light' ? 'Light' : 'Dark'}</div>
                      <div className={`text-xs mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                        {t === 'light' ? 'Clean and bright' : 'Easy on the eyes'}
                      </div>
                    </button>
                  ))}
                </div>
              </div>
            </div>

            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">Docker Configuration</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                  Container runtime settings
                </p>
              </div>
              <div className="p-6 space-y-4">
                <div>
                  <label className={labelClass}>GPU Acceleration</label>
                  <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-700' : 'bg-gray-50'}`}>
                    <label className="flex items-center gap-3 cursor-pointer">
                      <input
                        type="checkbox"
                        defaultChecked
                        className="w-5 h-5 text-blue-600 rounded border-gray-300 focus:ring-blue-500"
                      />
                      <div>
                        <p className="font-medium">Enable NVIDIA GPU</p>
                        <p className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                          Pass through NVIDIA GPU to containers for acceleration
                        </p>
                      </div>
                    </label>
                  </div>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* === EMAIL NOTIFICATIONS TAB === */}
        {activeTab === 'notifications' && (
          <div className="space-y-6">
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">Email Notifications</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                  Configure email notifications for job completion and alerts
                </p>
              </div>
              <div className="p-6 space-y-5">
                <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-700' : 'bg-gray-50'}`}>
                  <label className="flex items-center gap-3 cursor-pointer">
                    <input
                      type="checkbox"
                      checked={emailConfig.enabled}
                      onChange={(e) => setEmailConfig({ ...emailConfig, enabled: e.target.checked })}
                      className="w-5 h-5 text-blue-600 rounded"
                    />
                    <div>
                      <p className="font-medium">Enable Email Notifications</p>
                      <p className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                        Receive email alerts when jobs complete or fail
                      </p>
                    </div>
                  </label>
                </div>

                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <div>
                    <label className={labelClass}>SMTP Host</label>
                    <input
                      type="text"
                      value={emailConfig.smtpHost}
                      onChange={(e) => setEmailConfig({ ...emailConfig, smtpHost: e.target.value })}
                      className={inputClass()}
                      placeholder="smtp.gmail.com"
                    />
                  </div>
                  <div>
                    <label className={labelClass}>SMTP Port</label>
                    <input
                      type="number"
                      value={emailConfig.smtpPort}
                      onChange={(e) => setEmailConfig({ ...emailConfig, smtpPort: e.target.value })}
                      className={inputClass()}
                      placeholder="587"
                    />
                  </div>
                </div>

                <div>
                  <label className={labelClass}>From Email</label>
                  <input
                    type="email"
                    value={emailConfig.from}
                    onChange={(e) => setEmailConfig({ ...emailConfig, from: e.target.value })}
                    className={inputClass()}
                    placeholder="biodockify@example.com"
                  />
                </div>

                <div>
                  <label className={labelClass}>SMTP Password</label>
                  <input
                    type="password"
                    value={emailConfig.password}
                    onChange={(e) => setEmailConfig({ ...emailConfig, password: e.target.value })}
                    className={inputClass()}
                    placeholder="App-specific password"
                  />
                  <p className={hintClass}>Use an app-specific password for Gmail</p>
                </div>

                <div>
                  <label className={labelClass}>To Email</label>
                  <input
                    type="email"
                    value={emailConfig.to}
                    onChange={(e) => setEmailConfig({ ...emailConfig, to: e.target.value })}
                    className={inputClass()}
                    placeholder="your-email@example.com"
                  />
                </div>

                <div className="flex gap-3 pt-2">
                  <button
                    onClick={async () => {
                      try {
                        await fetch('/notifications/configure', {
                          method: 'POST',
                          headers: { 'Content-Type': 'application/json' },
                          body: JSON.stringify({
                            email_enabled: emailConfig.enabled,
                            email_smtp_host: emailConfig.smtpHost,
                            email_smtp_port: parseInt(emailConfig.smtpPort) || 587,
                            email_from: emailConfig.from,
                            email_password: emailConfig.password,
                            email_to: emailConfig.to,
                          }),
                        })
                        setMessage('Email settings saved')
                      } catch {
                        setMessage('Failed to save email settings')
                      }
                    }}
                    className="px-5 py-2.5 bg-blue-600 hover:bg-blue-700 text-white font-medium rounded-lg transition-colors"
                  >
                    Save Email Settings
                  </button>
                  <button
                    onClick={async () => {
                      try {
                        const res = await fetch('/notifications/test?channel=email', { method: 'POST' })
                        const data = await res.json()
                        setMessage(data.status === 'sent' ? 'Test email sent!' : data.message || 'Test failed')
                      } catch {
                        setMessage('Failed to send test email')
                      }
                    }}
                    className={`px-5 py-2.5 font-medium rounded-lg transition-colors ${
                      isDark ? 'bg-gray-700 hover:bg-gray-600 text-white border border-gray-600' : 'bg-white hover:bg-gray-50 text-gray-700 border border-gray-300'
                    }`}
                  >
                    Send Test Email
                  </button>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* === PLUGINS TAB === */}
        {activeTab === 'plugins' && (
          <div className="space-y-6">
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">Data Source Plugins</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                  Enable external data sources and services for CrewAI agents to fetch information from
                </p>
              </div>
              <div className="p-6 space-y-3">
                {plugins.map(plugin => (
                  <div
                    key={plugin.id}
                    className={`p-4 rounded-lg border ${
                      isDark ? 'bg-gray-700 border-gray-600' : 'bg-gray-50 border-gray-200'
                    }`}
                  >
                    <div className="flex items-center justify-between">
                      <div className="flex items-center gap-3">
                        <span className="text-2xl">{plugin.icon}</span>
                        <div>
                          <p className="font-medium">{plugin.name}</p>
                          <p className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                            {plugin.description}
                          </p>
                        </div>
                      </div>
                      <label className="relative inline-flex items-center cursor-pointer">
                        <input
                          type="checkbox"
                          checked={plugin.enabled}
                          onChange={() => {
                            setPlugins(prev => prev.map(p =>
                              p.id === plugin.id ? { ...p, enabled: !p.enabled } : p
                            ))
                          }}
                          className="sr-only peer"
                        />
                        <div className={`w-11 h-6 rounded-full peer peer-checked:after:translate-x-full after:content-[''] after:absolute after:top-[2px] after:left-[2px] after:bg-white after:rounded-full after:h-5 after:w-5 after:transition-all ${
                          isDark ? 'bg-gray-600 peer-checked:bg-blue-600' : 'bg-gray-300 peer-checked:bg-blue-600'
                        }`}>
                          <div className={`absolute top-[2px] left-[2px] bg-white border rounded-full h-5 w-5 transition-all ${
                            plugin.enabled ? 'translate-x-5' : ''
                          }`} />
                        </div>
                      </label>
                    </div>
                  </div>
                ))}
              </div>
            </div>
            
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">API Keys</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                  Configure API keys for external services
                </p>
              </div>
              <div className="p-6 space-y-4">
                <div>
                  <label className={labelClass}>PubChem API Key (optional)</label>
                  <input
                    type="password"
                    placeholder="Enter PubChem API key"
                    className={inputClass()}
                  />
                </div>
                <div>
                  <label className={labelClass}>ChEMBL API Key (optional)</label>
                  <input
                    type="password"
                    placeholder="Enter ChEMBL API key"
                    className={inputClass()}
                  />
                </div>
                <div>
                  <label className={labelClass}>ZINC Database Access (optional)</label>
                  <input
                    type="password"
                    placeholder="Enter ZINC access token"
                    className={inputClass()}
                  />
                </div>
              </div>
            </div>
            
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">Custom Plugins</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                  Add custom data source plugins
                </p>
              </div>
              <div className="p-6">
                <button className={`w-full py-3 rounded-lg border-2 border-dashed ${
                  isDark ? 'border-gray-600 text-gray-400 hover:border-blue-500 hover:text-blue-400' : 'border-gray-300 text-gray-500 hover:border-blue-500 hover:text-blue-600'
                }`}>
                  + Install Plugin
                </button>
              </div>
            </div>
          </div>
        )}

        {activeTab === 'about' && (
          <div className="space-y-6">
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">About BioDockify Studio AI</h2>
              </div>
              <div className="p-6">
                <div className="flex items-center gap-4 mb-6">
                  <div className="w-16 h-16 bg-gradient-to-br from-blue-500 to-purple-600 rounded-xl flex items-center justify-center text-2xl">
                    🧬
                  </div>
                  <div>
                    <h3 className="text-xl font-bold">BioDockify Studio AI</h3>
                    <p className={isDark ? 'text-gray-400' : 'text-gray-500'}>AI-Powered Drug Discovery Platform</p>
                  </div>
                </div>
                <div className="space-y-2">
                  <div className="flex justify-between py-2 border-b border-gray-200 dark:border-gray-700">
                    <span className={isDark ? 'text-gray-400' : 'text-gray-500'}>Version</span>
                    <span className="font-medium">{APP_VERSION}</span>
                  </div>
                  <div className="flex justify-between py-2 border-b border-gray-200 dark:border-gray-700">
                    <span className={isDark ? 'text-gray-400' : 'text-gray-500'}>Stack</span>
                    <span className="font-medium">React + FastAPI + RDKit + CrewAI</span>
                  </div>
                  <div className="flex justify-between py-2 border-b border-gray-200 dark:border-gray-700">
                    <span className={isDark ? 'text-gray-400' : 'text-gray-500'}>Docking</span>
                    <span className="font-medium">AutoDock Vina + GNINA + RF</span>
                  </div>
                  <div className="flex justify-between py-2">
                    <span className={isDark ? 'text-gray-400' : 'text-gray-500'}>License</span>
                    <span className="font-medium">MIT</span>
                  </div>
                </div>
              </div>
            </div>

            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">Features</h2>
              </div>
              <div className="p-6 grid grid-cols-2 gap-4">
                {[
                  { icon: '🧪', name: 'ChemDraw', desc: '2D/3D molecule editor' },
                  { icon: '🔬', name: 'Docking', desc: 'Vina/GNINA/RF scoring' },
                  { icon: '⚡', name: 'MD Simulation', desc: 'Molecular dynamics' },
                  { icon: '🤖', name: 'AI Assistant', desc: 'Multi-provider LLM chat' },
                  { icon: '📊', name: 'QSAR Modeling', desc: 'Predictive analysis' },
                  { icon: '💊', name: 'Pharmacophore', desc: 'Feature-based screening' },
                  { icon: '🛡️', name: 'ADMET', desc: 'Toxicity and drug-likeness' },
                  { icon: '📧', name: 'Notifications', desc: 'Email job alerts' },
                ].map(f => (
                  <div key={f.name} className={`p-3 rounded-lg ${isDark ? 'bg-gray-700' : 'bg-gray-50'}`}>
                    <div className="text-xl mb-1">{f.icon}</div>
                    <div className="font-medium text-sm">{f.name}</div>
                    <div className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>{f.desc}</div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}

        {activeTab === 'language' && (
          <div className="space-y-6">
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">🌐 Language</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Choose your preferred language</p>
              </div>
              <div className="p-6 grid grid-cols-2 gap-3">
                {locales.map(loc => (
                  <button
                    key={loc.code}
                    onClick={() => { setLocale(loc.code); setCurrentLocale(loc.code) }}
                    className={`p-4 rounded-lg border text-left transition-colors ${
                      currentLocale === loc.code
                        ? isDark ? 'bg-blue-900/50 border-blue-500' : 'bg-blue-50 border-blue-500'
                        : isDark ? 'bg-gray-700 border-gray-600 hover:border-gray-500' : 'bg-gray-50 border-gray-200 hover:border-gray-300'
                    }`}
                  >
                    <div className="font-medium">{loc.name}</div>
                    <div className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>{loc.code.toUpperCase()} • {loc.dir === 'rtl' ? 'RTL' : 'LTR'}</div>
                  </button>
                ))}
              </div>
            </div>
          </div>
        )}

        {activeTab === 'accessibility' && (
          <div className="space-y-6">
            <div className={`rounded-xl border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'} shadow-sm`}>
              <div className="p-6 border-b border-gray-200 dark:border-gray-700">
                <h2 className="text-lg font-semibold">♿ Accessibility</h2>
                <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Customize the interface for your needs</p>
              </div>
              <div className="p-6 space-y-4">
                <label className="flex items-center gap-3 cursor-pointer">
                  <input type="checkbox" checked={highContrast} onChange={toggleHighContrast} className="w-5 h-5 rounded" />
                  <div>
                    <p className="font-medium">High Contrast Mode</p>
                    <p className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Increase contrast for better visibility</p>
                  </div>
                </label>
                <label className="flex items-center gap-3 cursor-pointer">
                  <input type="checkbox" checked={reducedMotion} onChange={toggleReducedMotion} className="w-5 h-5 rounded" />
                  <div>
                    <p className="font-medium">Reduced Motion</p>
                    <p className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Minimize animations and transitions</p>
                  </div>
                </label>
                <div>
                  <p className="font-medium mb-2">Font Size: {fontSize}px</p>
                  <input type="range" min="12" max="24" value={fontSize} onChange={e => setFontSize(Number(e.target.value))} className="w-full" />
                  <div className="flex justify-between text-xs text-gray-500 mt-1">
                    <span>Small (12px)</span>
                    <span>Large (24px)</span>
                  </div>
                </div>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}
