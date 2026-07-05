#!/usr/bin/env python3
"""
eval_plant_virus — 统一基准测试入口
=====================================
用法:
  python run_benchmark.py status              # 查看各阶段完成状态
  python run_benchmark.py eval1               # 运行评估一（宿主过滤）
  python run_benchmark.py eval2               # 运行评估二（已知病毒检测）
  python run_benchmark.py all                 # 按依赖顺序运行全部
  python run_benchmark.py all --auto-prereq   # 自动运行缺失的前置步骤

阶段依赖:
  prep ─┬─→ eval1 (宿主过滤)
        ├─→ eval2 (已知病毒检测)
        ├─→ eval3 (病毒组装)
        ├─→ eval4 (候选病毒鉴定, 需 prep_master_eval)
        ├─→ eval5 (病毒分类, 需 prep_build_class_eval)
        ├─→ eval6 (序列去重聚类)
        └─→ eval7 (宿主分类, 需 eval5 输出)
"""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# Windows GBK 编码兼容
try:
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass  # Python < 3.7 或 stdin 重定向时不支持

PROJECT_ROOT = Path(__file__).resolve().parent
EVAL_DIR = PROJECT_ROOT / "eval"

# ── 阶段定义: (shell脚本, 前置依赖列表, 输出目录, DONE标记) ──
PHASES = {
    "prep": {
        "desc": "选取评估病毒 (Step 1)",
        "script": None,  # 手动运行
        "prereqs": [],
        "done_marker": "step1_eval_viruses/.DONE",
        "output_dir": "step1_eval_viruses",
    },
    "simulate": {
        "desc": "模拟测序数据 (Step 2-4)",
        "script": None,  # 手动运行
        "prereqs": ["prep"],
        "done_marker": "step2_simulator/.DONE",
        "output_dir": "step2_simulator",
    },
    "eval1": {
        "desc": "宿主过滤消融",
        "script": "run_eval_host_depletion.sh",
        "prereqs": ["simulate"],
        "done_marker": "step5_host_free/.DONE",
        "output_dir": "step5_host_free",
    },
    "eval2": {
        "desc": "已知病毒检测",
        "script": "run_eval_known_virus.sh",
        "prereqs": ["simulate", "prep"],
        "done_marker": "step6_identify_PVDB/.DONE",
        "output_dir": "step6_identify_PVDB",
    },
    "eval3": {
        "desc": "病毒组装",
        "script": "run_eval_assembly.sh",
        "prereqs": ["simulate", "prep"],
        "done_marker": "step7_assembly/.DONE",
        "output_dir": "step7_assembly",
    },
    "eval4": {
        "desc": "候选病毒鉴定",
        "script": "run_eval_identification.sh",
        "prereqs": ["prep"],
        "done_marker": "step8_result/.DONE",
        "output_dir": "step8_result",
    },
    "eval5": {
        "desc": "病毒分类",
        "script": "run_eval_classification.sh",
        "prereqs": ["prep"],
        "done_marker": "step9_classification/.DONE",
        "output_dir": "step9_classification",
    },
    "eval6": {
        "desc": "序列去重聚类",
        "script": "run_dedup_clustering.sh",
        "prereqs": ["prep"],
        "done_marker": "step10_dedup/.DONE",
        "output_dir": "step10_dedup",
    },
    "eval7": {
        "desc": "宿主分类基准",
        "script": "run_eval_host_prediction.sh",
        "prereqs": ["eval5"],
        "done_marker": "step11_host_evaluation/.DONE",
        "output_dir": "step11_host_evaluation",
    },
}

STATE_FILE = PROJECT_ROOT / ".benchmark_state.json"


# ── 工具函数 ────────────────────────────────────────────────

def load_state():
    if STATE_FILE.exists():
        return json.loads(STATE_FILE.read_text())
    return {}


def save_state(state):
    STATE_FILE.write_text(json.dumps(state, indent=2))


def check_phase(phase_name):
    """检查阶段是否完成"""
    info = PHASES[phase_name]
    marker = PROJECT_ROOT / info["done_marker"]
    return marker.exists()


def check_prereqs(phase_name):
    """检查前置依赖是否满足，返回缺失列表"""
    info = PHASES[phase_name]
    return [p for p in info["prereqs"] if not check_phase(p)]


def run_phase(phase_name, auto_prereq=False, extra_args=None):
    """运行单个评估阶段"""
    info = PHASES[phase_name]

    if check_phase(phase_name):
        print(f"[{phase_name}] ✅ 已完成 ({info['desc']})")
        return True

    # 检查前置依赖
    missing = check_prereqs(phase_name)
    if missing:
        if auto_prereq:
            print(f"[{phase_name}] ⚠️  前置阶段未完成: {missing}，自动运行...")
            for p in missing:
                if not run_phase(p, auto_prereq=True):
                    print(f"[{phase_name}] ❌ 前置 {p} 失败")
                    return False
        else:
            print(f"[{phase_name}] ❌ 前置阶段未完成: {missing}")
            print(f"   请先运行: python run_benchmark.py {' '.join(missing)}")
            print(f"   或使用 --auto-prereq 自动运行")
            return False

    if info["script"] is None:
        print(f"[{phase_name}] ⚠️  无自动脚本，请手动运行: {info['desc']}")
        return True

    script_path = EVAL_DIR / info["script"]
    if not script_path.exists():
        print(f"[{phase_name}] ❌ 脚本不存在: {script_path}")
        return False

    # 运行
    print(f"\n{'='*60}")
    print(f"[{phase_name}] 🚀 {info['desc']}")
    print(f"[{phase_name}] 脚本: {script_path}")
    start = time.time()

    cmd = ["bash", str(script_path)]
    if extra_args:
        cmd.extend(extra_args)

    try:
        result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
        elapsed = time.time() - start

        if result.returncode == 0:
            print(f"[{phase_name}] ✅ 完成 ({elapsed:.0f}s)")
            # 更新状态
            state = load_state()
            state[phase_name] = {
                "completed": datetime.now().isoformat(),
                "duration_s": round(elapsed, 1),
                "output_dir": info["output_dir"],
            }
            save_state(state)
            return True
        else:
            print(f"[{phase_name}] ❌ 失败 (exit={result.returncode}, {elapsed:.0f}s)")
            return False
    except KeyboardInterrupt:
        print(f"\n[{phase_name}] ⏹️  用户中断")
        return False


# ── 命令实现 ────────────────────────────────────────────────

def cmd_status(args):
    """显示各阶段状态"""
    print("\n评估阶段状态:")
    print(f"{'阶段':<10} {'状态':<10} {'描述'}")
    print("-" * 60)
    state = load_state()
    for name, info in PHASES.items():
        done = check_phase(name)
        status_str = "✅ 完成" if done else "⬜ 待运行"
        print(f"{name:<10} {status_str:<10} {info['desc']}")
    print("-" * 60)
    print(f"状态文件: {STATE_FILE}")


def cmd_run(args):
    """运行指定阶段"""
    phases_to_run = args.phases
    if "all" in phases_to_run:
        phases_to_run = list(PHASES.keys())

    for phase in phases_to_run:
        if phase not in PHASES:
            print(f"❌ 未知阶段: {phase}")
            print(f"   有效阶段: {', '.join(PHASES.keys())}")
            sys.exit(1)

    for phase in phases_to_run:
        if not run_phase(phase, auto_prereq=args.auto_prereq, extra_args=args.args):
            if not args.keep_going:
                print(f"\n❌ {phase} 失败，终止。使用 --keep-going 忽略错误继续。")
                sys.exit(1)


def cmd_validate(args):
    """验证输出"""
    phase = args.phase if args.phase != "all" else None
    phases_to_check = [phase] if phase else list(PHASES.keys())

    all_ok = True
    for p in phases_to_check:
        info = PHASES[p]
        output_dir = PROJECT_ROOT / info["output_dir"]
        done_marker = PROJECT_ROOT / info["done_marker"]

        if done_marker.exists():
            print(f"  ✅ {p:<8} DONE marker: {info['done_marker']}")
        else:
            print(f"  ❌ {p:<8} 缺少 DONE: {info['done_marker']}")
            all_ok = False

        if output_dir.exists():
            file_count = sum(1 for _ in output_dir.rglob("*") if _.is_file())
            print(f"     {p:<8} 输出目录: {info['output_dir']} ({file_count} files)")
        else:
            print(f"     {p:<8} 输出目录不存在: {info['output_dir']}")
            all_ok = False

    sys.exit(0 if all_ok else 1)


# ── CLI ─────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="eval_plant_virus 统一基准测试入口",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python run_benchmark.py status              # 查看状态
  python run_benchmark.py eval1               # 运行评估一
  python run_benchmark.py eval1 eval2         # 运行多个
  python run_benchmark.py all --auto-prereq   # 全自动运行
  python run_benchmark.py validate            # 验证输出
        """,
    )

    sub = parser.add_subparsers(dest="command", help="命令")

    # status
    p_status = sub.add_parser("status", help="查看各阶段完成状态")

    # run
    p_run = sub.add_parser("run", help="运行评估阶段")
    p_run.add_argument("phases", nargs="+",
                       choices=list(PHASES.keys()) + ["all"],
                       help="要运行的阶段 (all = 全部)")
    p_run.add_argument("--auto-prereq", action="store_true",
                       help="自动运行缺失的前置阶段")
    p_run.add_argument("--keep-going", action="store_true",
                       help="某个阶段失败后继续运行后续")
    p_run.add_argument("args", nargs="*", help="传递给 shell 脚本的额外参数")

    # validate
    p_val = sub.add_parser("validate", help="验证输出完整性")
    p_val.add_argument("phase", nargs="?", default="all",
                       choices=list(PHASES.keys()) + ["all"],
                       help="要验证的阶段 (默认: all)")

    # 兼容: 直接传阶段名 (如 python run_benchmark.py eval1)
    if len(sys.argv) >= 2 and sys.argv[1] in list(PHASES.keys()) + ["all"]:
        # 重新构造参数
        phases = [sys.argv[1]]
        extra = sys.argv[2:]
        sys.argv = [sys.argv[0], "run"] + phases + extra

    args = parser.parse_args()

    if args.command == "status":
        cmd_status(args)
    elif args.command == "run":
        cmd_run(args)
    elif args.command == "validate":
        cmd_validate(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
