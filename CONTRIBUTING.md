# Contributing to BioDockify Docking Studio

Thank you for your interest in contributing to BioDockify Docking Studio!

## Code of Conduct

By participating in this project, you agree to maintain a respectful and inclusive environment for everyone.

## Ways to Contribute

### Reporting Bugs

1. Check if the bug has already been reported
2. Use the issue template to create a bug report
3. Include steps to reproduce, expected behavior, and actual behavior
4. Include your environment details (OS, Python version, Docker version)

### Suggesting Features

1. Check if the feature has been discussed before
2. Describe the feature in detail
3. Explain why this feature would be useful
4. Provide use cases

### Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests if applicable
5. Ensure all tests pass
6. Commit your changes with clear commit messages
7. Push to your fork
8. Submit a pull request

## Development Setup

### Prerequisites

- Python 3.9+
- Docker & Docker Compose
- Git

### Local Development

```bash
# Clone the repository
git clone https://github.com/tajo9128/Docking-studio.git
cd Docking-studio

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or
venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Start Docker services
docker-compose up -d

# Run the application
python -m src.biodockify_main
```

### Running Tests

```bash
# Run Python tests
pytest tests/

# Run security scans
docker-compose run --rm backend bandit -r src/
docker-compose run --rm backend safety check
```

## Coding Standards

### Python

- Follow PEP 8
- Use type hints where possible
- Add docstrings to functions and classes
- Keep functions small and focused

### Git Commit Messages

- Use clear, descriptive commit messages
- Start with a verb (Add, Fix, Update, Remove)
- Reference issues where applicable

Example:
```
Add SSE streaming for real-time docking progress

- Implement /dock/stream endpoint
- Add DockingWorker class
- Update progress bar UI
Fixes #123
```

### UI/UX Guidelines

- Follow the existing design system in `theme.py`
- Use Design Tokens for colors, spacing, typography
- Test on multiple screen sizes

## Documentation

- Update README.md for user-facing changes
- Add docstrings to new functions
- Update API documentation if endpoints change

## Recognition

Contributors will be recognized in:
- README.md contributors section
- CHANGELOG.md
- Release notes

## License

By contributing, you agree that your contributions will be licensed under the Apache License 2.0.
