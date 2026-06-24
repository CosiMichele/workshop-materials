#!/bin/bash

# Add Claude API environment variables to bashrc

cat >> ~/.bashrc << 'EOF'

# Claude API Configuration
export ANTHROPIC_AUTH_TOKEN="insert-key-here"
export ANTHROPIC_BEDROCK_BASE_URL=https://llm-api.cyverse.ai/bedrock
export CLAUDE_CODE_SKIP_BEDROCK_AUTH=1
export CLAUDE_CODE_USE_BEDROCK=1
export ANTHROPIC_DEFAULT_OPUS_MODEL="ai2s-claude-opus-4-6"
export ANTHROPIC_DEFAULT_SONNET_MODEL="ai2s-claude-sonnet-4-6"
export ANTHROPIC_DEFAULT_HAIKU_MODEL="ai2s-claude-haiku-4-5"

EOF

# Source bashrc to apply changes immediately
source ~/.bashrc

echo "? Environment variables added to ~/.bashrc"
echo "? Changes applied to current session"
