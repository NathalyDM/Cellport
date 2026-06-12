from pathlib import Path

import json

base = Path('.')
notebooks = list(base.glob('*.ipynb')) + list(base.glob('Cargo Segmentation/*.ipynb'))
replacements = {
    "C:/Users/utraf": "C:/Users/<anon>",
    "C:\\Users\\utraf": "C:\\Users\\<anon>",
    "Clinical_AI": "<anon>",
}

def sanitize_notebook(path: Path):
    try:
        nb = json.loads(path.read_text(encoding='utf-8'))
    except Exception:
        print('Skipping (not valid JSON):', path)
        return
    # Clear outputs and execution counts
    for cell in nb.get('cells', []):
        if isinstance(cell, dict):
            cell['outputs'] = []
            cell['execution_count'] = None
            # sanitize source lines
            if 'source' in cell:
                cell['source'] = [s.replace(k, v) for s in cell['source'] for k, v in replacements.items()]
    # Minimal metadata
    nb['metadata'] = {}
    content = json.dumps(nb, ensure_ascii=False, indent=1)
    for old, new in replacements.items():
        content = content.replace(old, new)
    path.write_text(content, encoding='utf-8')
    print('Sanitized:', path)

for nb_path in notebooks:
    sanitize_notebook(Path(nb_path))

print('Done.')
