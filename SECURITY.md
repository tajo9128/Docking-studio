# Security Policy

## Supported Versions

We release security updates for the following versions:

| Version | Supported          |
| ------- | ------------------ |
| 1.2.x   | :white_check_mark: |
| 1.1.x   | :white_check_mark: |
| 1.0.x   | :x:                |

## Reporting a Vulnerability

If you discover a security vulnerability within BioDockify Docking Studio, please send an email to the development team. All security vulnerabilities will be promptly addressed.

### What to Include

- Type of vulnerability
- Full paths of source file(s) related to the vulnerability
- Location of the affected source code
- Any special configuration required to reproduce the issue
- Step-by-step instructions to reproduce the issue
- Proof-of-concept or exploit code (if possible)
- Impact of the issue

## Security Features

### Built-in Security

- **Input Validation**: All user inputs are validated and sanitized
- **Secure File Handling**: Files are validated before processing
- **SQL Injection Protection**: Parameterized queries for database operations
- **CORS Protection**: Configured to restrict cross-origin requests
- **Rate Limiting**: Prevents abuse of API endpoints

### Docker Security

- Non-root user execution
- Minimal base image
- Read-only filesystem where possible
- Regular security scanning with Trivy

### AI Safety

- Ollama runs in isolated Docker container
- No external model downloads without user consent
- Offline mode available for air-gapped environments

## Dependencies Security

We regularly scan our dependencies for vulnerabilities:

- **Trivy**: Container vulnerability scanning
- **Bandit**: Python security linting
- **Safety**: Python dependency checking
- **Gitleaks**: Secrets detection

## Compliance

For enterprise deployments, the system supports:

- Air-gapped operation (offline mode)
- Local model execution
- No cloud dependencies required
- Complete audit logging

## Contact

For security issues, please do not open a public GitHub issue. Contact the maintainers directly through the repository.
