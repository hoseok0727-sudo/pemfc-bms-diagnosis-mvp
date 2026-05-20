from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "ui_recordings"


def render_dashboard_gif(
    dataset: str,
    output_path: str | Path,
    base_url: str = "http://localhost:8765/dashboard.html",
    frames: int = 28,
    duration_ms: int = 130,
    width: int = 1280,
    height: int = 820,
) -> Path:
    chrome = _find_chrome()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="pemfc_dashboard_frames_") as temp_dir:
        temp_dir_path = Path(temp_dir)
        frame_paths: list[Path] = []
        for i in range(frames):
            pct = i / (frames - 1)
            url = f"{base_url}?dataset={dataset}&frame={pct:.5f}&autoplay=0"
            png_path = temp_dir_path / f"{dataset}_{i:03d}.png"
            cmd = [
                str(chrome),
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--no-first-run",
                f"--window-size={width},{height}",
                "--hide-scrollbars",
                "--run-all-compositor-stages-before-draw",
                "--virtual-time-budget=2500",
                f"--screenshot={png_path}",
                url,
            ]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            frame_paths.append(png_path)

        images = [Image.open(path).convert("P", palette=Image.Palette.ADAPTIVE) for path in frame_paths]
        images[0].save(
            output_path,
            save_all=True,
            append_images=images[1:],
            optimize=True,
            duration=duration_ms,
            loop=0,
        )
        for image in images:
            image.close()
    return output_path


def render_default_gifs(output_dir: str | Path = DEFAULT_OUTPUT_DIR) -> list[Path]:
    output_dir = Path(output_dir)
    return [
        render_dashboard_gif("normal_data", output_dir / "normal_operation_dashboard.gif"),
        render_dashboard_gif("bad_data", output_dir / "fault_critical_dashboard.gif"),
    ]


def _find_chrome() -> Path:
    candidates = [
        Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"),
        Path(r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"),
        Path(r"C:\Program Files\Microsoft\Edge\Application\msedge.exe"),
        Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    resolved = shutil.which("chrome") or shutil.which("msedge")
    if resolved:
        return Path(resolved)
    raise FileNotFoundError("Chrome or Edge executable was not found.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render dashboard playback GIFs.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--frames", type=int, default=28)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    normal = render_dashboard_gif(
        "normal_data",
        args.output_dir / "normal_operation_dashboard.gif",
        frames=args.frames,
    )
    fault = render_dashboard_gif(
        "bad_data",
        args.output_dir / "fault_critical_dashboard.gif",
        frames=args.frames,
    )
    print(f"Normal operation GIF saved to: {normal}")
    print(f"Fault/Critical GIF saved to: {fault}")


if __name__ == "__main__":
    main()
