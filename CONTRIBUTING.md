# Contributing to VecMap

We welcome contributions to VecMap. This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone https://github.com/YOUR_USERNAME/VecMap.git`
3. Create a new branch: `git checkout -b feature-name`
4. Make your changes
5. Run tests: `pytest`
6. Commit with clear messages: `git commit -m "Add feature X"`
7. Push to your fork: `git push origin feature-name`
8. Create a pull request

## Development Setup

```bash
git clone https://github.com/the-jordan-lab/VecMap.git
cd VecMap
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e ".[dev]"
```

## Code Style

- Follow PEP 8
- Use type hints where appropriate
- Add docstrings to all public functions
- Keep the core algorithm simple and readable

## Testing

- Add tests for new features
- Ensure all tests pass before submitting PR
- Include benchmark results for performance-related changes

## Pull Request Process

1. Update the README.md with details of changes if needed
2. Update the version number in `vecmap/__init__.py` if appropriate
3. The PR will be merged once reviewed and approved

## Code of Conduct

Please be respectful and professional in all interactions. 