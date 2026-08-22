"""
Workspace Builder Module for ATobb Combine Analysis.
Compiles Combine datacards (.txt) into RooWorkspaces (.root) using text2workspace.py.
"""

import os
import subprocess


def build_workspace(
    datacard_path,
    workspace_path=None,
    work_dir=None,
    physics_model="modules.charged_higgs_br:chargedHiggsBR",
    channel_masks=True,
):
    """
    Compile a datacard text file into a RooWorkspace ROOT file using text2workspace.py.
    """
    if not os.path.exists(datacard_path):
        print(f"Error: Datacard file {datacard_path} not found.")
        return False

    if workspace_path is None:
        workspace_path = os.path.splitext(datacard_path)[0] + "_ws.root"

    ws_dir = os.path.dirname(workspace_path)
    if ws_dir:
        os.makedirs(ws_dir, exist_ok=True)

    env = os.environ.copy()
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cur_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{repo_root}:{cur_pythonpath}" if cur_pythonpath else repo_root

    model_flag = f"-P {physics_model}" if physics_model else ""
    mask_flag = "--channel-masks" if channel_masks else ""
    cmd = f"text2workspace.py {datacard_path} -o {workspace_path} {model_flag} {mask_flag}"

    res = subprocess.run(
        cmd,
        shell=True,
        cwd=work_dir,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )

    if res.returncode != 0:
        print(f"[-] text2workspace.py failed for {datacard_path}:\n{res.stderr}")
        return False

    return workspace_path
