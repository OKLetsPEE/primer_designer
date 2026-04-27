import re
import unicodedata

import primer3
import streamlit as st


# --- 页面配置 ---
st.set_page_config(page_title="Chanz觉得你的时间很宝贵")

st.title("Chanz帮你设计引物")
st.markdown(
    "粘贴序列前请先在左侧边栏选择设计模式 | "
    "常规引物设计请选择模式3 | 模式 1 & 2 为 circRNA 专用 | Written by Chanz"
)


# =========================
# 输入清洗与校验函数
# =========================
def _longest_base_run(text: str) -> int:
    """计算未被空格/数字/分隔符打断的最长 A/C/G/T/U/N 连续长度。"""
    runs = re.findall(r"[ACGTUNacgtun]+", text)
    return max((len(run) for run in runs), default=0)


def extract_nucleotide_sequence(raw_text: str, min_candidate_len: int = 10):
    """
    从用户输入的混杂文本中提取核酸序列。

    设计目标：
    1. 允许用户粘贴 FASTA、带行号的序列、带空格/换行/数字的序列。
    2. 自动忽略说明性文字，例如：xxx mRNA sequenced、gene name、description 等。
    3. 避免把英文说明文字中的零散 A/C/G/T 误当作真实序列。
    4. 允许 U，并自动转换为 T；允许 N，但会提示用户注意。

    核心思路：逐行识别“像序列的片段”，而不是简单保留所有 A/C/G/T。
    这样可以避免把 mRNA、sequence、annotation 等英文单词里的碱基字母误并入模板。

    返回：clean_seq, report
    """
    report = {
        "raw_len": len(raw_text or ""),
        "headers_skipped": 0,
        "candidate_blocks": 0,
        "short_blocks_ignored": 0,
        "u_to_t": 0,
        "n_count": 0,
        "fallback_used": False,
        "warnings": [],
    }

    if not raw_text:
        return "", report

    # 统一全角/半角字符，处理从网页、Word、PDF 复制来的奇怪字符
    text = unicodedata.normalize("NFKC", raw_text)
    text = text.replace("\ufeff", "")

    # 只把“由 A/C/G/T/U/N + 数字 + 空白 + 常见分隔符”组成的连续片段视为候选序列块。
    # 英文说明中的其他字母会自然切断候选块，从而减少误提取。
    # 注意：这里不跨行匹配，避免上一行说明文字中的 NA/RNA 与下一行真实序列被错误拼接。
    candidate_pattern = re.compile(r"[ACGTUNacgtun\d \t\-_.|/\\]+")

    candidates = []

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue

        # 跳过 FASTA 标题行，避免标题里的 A/C/G/T 被误并入序列
        if stripped.startswith(">"):
            report["headers_skipped"] += 1
            continue

        # 如果一行中存在 “Sequence: ATGC...” / “seq = ATGC...” 这类格式，优先分析分隔符之后的内容。
        # 同时保留整行作为备选，以兼容没有冒号/等号的行。
        segments = [line]
        if ":" in line:
            segments.insert(0, line.split(":")[-1])
        if "=" in line:
            segments.insert(0, line.split("=")[-1])

        accepted_this_line = False
        for segment in segments:
            letters_in_segment = re.findall(r"[A-Za-z]", segment)
            for match in candidate_pattern.finditer(segment):
                block = match.group(0)
                seq = re.sub(r"[^ACGTUNacgtun]", "", block).upper()
                if not seq:
                    continue

                # 判定这个片段是否真的“像序列”：
                # 1) 有足够长的连续碱基串；或
                # 2) 这一段文本里绝大多数英文字母都是碱基，适合处理 A T G C 或 ATG CCG 这种带空格格式。
                longest_run = _longest_base_run(block)
                base_ratio = len(seq) / len(letters_in_segment) if letters_in_segment else 0
                block_density = len(seq) / len(block.strip()) if block.strip() else 0
                looks_like_sequence = (
                    longest_run >= min_candidate_len
                    or base_ratio >= 0.70
                    or (len(seq) >= 20 and block_density >= 0.55)
                )

                # 很多 FASTA/网页复制序列的最后一行会短于 10 bp。
                # 如果前面已经识别到序列，且当前短片段本身很像序列，则把它作为尾部片段保留。
                short_tail_of_existing_sequence = (
                    len(seq) < min_candidate_len
                    and bool(candidates)
                    and base_ratio >= 0.70
                    and longest_run >= 4
                )

                if not looks_like_sequence and not short_tail_of_existing_sequence:
                    report["short_blocks_ignored"] += 1
                    continue

                if len(seq) < min_candidate_len and not short_tail_of_existing_sequence:
                    report["short_blocks_ignored"] += 1
                    continue

                report["candidate_blocks"] += 1
                report["u_to_t"] += seq.count("U")
                candidates.append(seq.replace("U", "T"))
                accepted_this_line = True

            # 如果分隔符后的 segment 已经提取成功，就不再重复分析整行，避免重复加入同一段序列
            if accepted_this_line:
                break

    clean_seq = "".join(candidates)

    # 兜底逻辑：有些用户可能输入的是非常碎片化的序列，例如 A T G C T A ...
    # 如果主逻辑没有提取到序列，但整段文本中绝大多数英文字母都是碱基，则启用宽松模式。
    if not clean_seq:
        # 去掉 FASTA 标题再做兜底，仍然避免标题污染
        no_header_lines = [line for line in text.splitlines() if not line.strip().startswith(">")]
        text_no_header = "\n".join(no_header_lines)
        letters = re.findall(r"[A-Za-z]", text_no_header)
        loose_seq = "".join(re.findall(r"[ACGTUacgtu]", text_no_header)).upper()
        base_ratio = len(loose_seq) / len(letters) if letters else 0
        if len(loose_seq) >= min_candidate_len and base_ratio >= 0.85:
            report["fallback_used"] = True
            report["u_to_t"] += loose_seq.count("U")
            clean_seq = loose_seq.replace("U", "T")

    report["n_count"] = clean_seq.count("N")

    if report["headers_skipped"]:
        report["warnings"].append(f"已自动忽略 {report['headers_skipped']} 行 FASTA 标题。")
    if report["u_to_t"]:
        report["warnings"].append(f"检测到 {report['u_to_t']} 个 U，已按 DNA 模板自动转换为 T。")
    if report["n_count"]:
        report["warnings"].append(
            f"清洗后的序列中仍有 {report['n_count']} 个 N；Primer3 通常会避开 N 区域，但 N 太多可能导致无结果。"
        )
    if report["short_blocks_ignored"]:
        report["warnings"].append(
            "已忽略若干过短或不像真实序列的片段；这通常来自说明性英文中的零散 A/C/G/T。"
        )
    if report["fallback_used"]:
        report["warnings"].append("已启用宽松清洗模式：输入整体看起来像核酸序列，因此按 A/T/G/C/U 提取。")

    return clean_seq, report


def calculate_gc(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def build_primer3_args(clean_seq: str, mode: int, primer_settings: dict | None = None):
    """根据模式生成 primer3 的 sequence args 与 global args。"""
    seq_len = len(clean_seq)

    global_args = {
        "PRIMER_NUM_RETURN": 5,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 57.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_MIN_GC": 40.0,
        "PRIMER_MAX_GC": 60.0,
        "PRIMER_MAX_POLY_X": 5,
    }

    if primer_settings:
        global_args.update(primer_settings)

    final_sequence = clean_seq
    seq_args = {
        "SEQUENCE_ID": "Target",
        "SEQUENCE_TEMPLATE": final_sequence,
    }
    status_msg = ""

    if mode == 1:
        if seq_len < 20:
            raise ValueError("模式1需要至少 20 bp 以上的输入序列，否则无法可靠设置 BSJ 附近的靶区。")
        final_sequence = clean_seq + clean_seq
        join_point = seq_len
        seq_args["SEQUENCE_TEMPLATE"] = final_sequence
        # SEQUENCE_TARGET 属于 sequence args，不建议放在 global args 中
        seq_args["SEQUENCE_TARGET"] = [join_point - 10, 20]
        global_args["PRIMER_PRODUCT_SIZE_RANGE"] = [[80, 150]]
        status_msg = "模式1：已执行序列串联模拟环化，并强制产物跨越 BSJ 附近区域。"

    elif mode == 2:
        if seq_len < 20:
            raise ValueError("模式2需要至少 20 bp 以上的输入序列，否则无法可靠设置 BSJ 附近的靶区。")
        final_sequence = clean_seq + clean_seq
        join_point = seq_len
        seq_args["SEQUENCE_TEMPLATE"] = final_sequence
        seq_args["SEQUENCE_TARGET"] = [join_point - 5, 10]
        global_args["PRIMER_PRODUCT_SIZE_RANGE"] = [[250, 600], [200, 800]]
        status_msg = "模式2：Sanger 测序模式，已执行序列串联模拟环化，并强制产物跨越 BSJ 附近区域。"

    elif mode == 3:
        global_args["PRIMER_PRODUCT_SIZE_RANGE"] = [[80, 200]]
        status_msg = "模式3：常规线性引物设计。"

    return seq_args, global_args, status_msg


# --- 侧边栏：参数设置 ---
with st.sidebar:
    st.header("参数设置")
    mode = st.radio(
        "选择设计模式",
        (1, 2, 3),
        format_func=lambda x: {
            1: "1. circRNA Divergent Primer (For qPCR)",
            2: "2. circRNA BSJ Primer (For Sanger-seq)",
            3: "3. 常规线性 mRNA Primer (For qPCR)",
        }[x],
    )
    st.info(
        "模式说明：\n\n"
        "**模式1**: 模拟环化，强制跨接头，产物 80-150bp\n\n"
        "**模式2**: 模拟环化，强制跨接头，产物 250-600bp\n\n"
        "**模式3**: 常规线性设计，产物 80-200bp"
    )

    with st.expander("高级参数：Tm / GC / 引物长度", expanded=False):
        st.caption(
            "默认参数会强烈优先选择 Tm 接近 60°C 的引物；"
            "如果你想看到更多变化，可以在这里调整。"
        )
        primer_opt_tm = st.number_input("目标 Tm (°C)", min_value=45.0, max_value=75.0, value=60.0, step=0.5)
        tm_col1, tm_col2 = st.columns(2)
        with tm_col1:
            primer_min_tm = st.number_input("最低 Tm (°C)", min_value=40.0, max_value=75.0, value=57.0, step=0.5)
        with tm_col2:
            primer_max_tm = st.number_input("最高 Tm (°C)", min_value=40.0, max_value=80.0, value=63.0, step=0.5)

        size_col1, size_col2, size_col3 = st.columns(3)
        with size_col1:
            primer_min_size = st.number_input("最短引物 nt", min_value=15, max_value=35, value=18, step=1)
        with size_col2:
            primer_opt_size = st.number_input("最佳引物 nt", min_value=15, max_value=35, value=20, step=1)
        with size_col3:
            primer_max_size = st.number_input("最长引物 nt", min_value=15, max_value=35, value=25, step=1)

        gc_col1, gc_col2 = st.columns(2)
        with gc_col1:
            primer_min_gc = st.number_input("最低 GC (%)", min_value=20.0, max_value=80.0, value=40.0, step=1.0)
        with gc_col2:
            primer_max_gc = st.number_input("最高 GC (%)", min_value=20.0, max_value=90.0, value=60.0, step=1.0)

    primer_settings = {
        "PRIMER_OPT_TM": float(primer_opt_tm),
        "PRIMER_MIN_TM": float(primer_min_tm),
        "PRIMER_MAX_TM": float(primer_max_tm),
        "PRIMER_MIN_SIZE": int(primer_min_size),
        "PRIMER_OPT_SIZE": int(primer_opt_size),
        "PRIMER_MAX_SIZE": int(primer_max_size),
        "PRIMER_MIN_GC": float(primer_min_gc),
        "PRIMER_MAX_GC": float(primer_max_gc),
    }

# --- 主界面 ---
st.subheader("输入序列")
raw_input = st.text_area(
    "在此粘贴基因序列 (5' → 3')",
    height=180,
    placeholder=(
        ">GeneA mRNA sequence\n"
        "1 ATG CCG TTA GCT ... 60\n"
        "也可以包含说明文字、数字、空格和换行，程序会自动清洗。"
    ),
)

uploaded_file = st.file_uploader(
    "也可以上传包含序列的文本文件 / FASTA 文件 / CSV 文件",
    type=["txt", "fa", "fasta", "fna", "csv"],
)

file_text = ""
if uploaded_file is not None:
    try:
        file_text = uploaded_file.read().decode("utf-8-sig")
    except UnicodeDecodeError:
        uploaded_file.seek(0)
        file_text = uploaded_file.read().decode("gbk", errors="ignore")

# 如果既粘贴又上传，则合并处理，方便用户在上传文件后临时补充序列片段
combined_input = "\n".join(part for part in [raw_input, file_text] if part)
clean_seq, clean_report = extract_nucleotide_sequence(combined_input)

if combined_input:
    st.caption(
        f"清洗后可用于设计的序列长度：{len(clean_seq)} bp "
        f"| GC：{calculate_gc(clean_seq):.1f}%"
    )
    if clean_report["warnings"]:
        for msg in clean_report["warnings"]:
            st.info(msg)

    with st.expander("查看清洗后的序列", expanded=False):
        if clean_seq:
            st.code(clean_seq, language="text")
            st.download_button(
                "下载清洗后的序列 FASTA",
                data=f">cleaned_sequence\n{clean_seq}\n",
                file_name="cleaned_sequence.fasta",
                mime="text/plain",
            )
        else:
            st.warning("暂未从输入中识别到足够长的核酸序列片段。")


if st.button("开始设计引物"):
    if not clean_seq:
        st.error("没有识别到有效核酸序列。请粘贴或上传包含 A/T/G/C 的序列后再试。")
        st.stop()

    if len(clean_seq) < 80:
        st.warning(
            "当前清洗后的序列长度低于 80 bp。由于产物范围通常从 80 bp 起，"
            "Primer3 可能找不到合适引物。"
        )

    if clean_report["n_count"] / len(clean_seq) > 0.05:
        st.warning("N 的比例超过 5%，建议尽量使用更完整的 A/T/G/C 序列，否则可能影响引物设计结果。")

    try:
        seq_args, global_args, status_msg = build_primer3_args(clean_seq, mode, primer_settings)
        results = primer3.bindings.design_primers(seq_args, global_args)
        pairs_found = results.get("PRIMER_PAIR_NUM_RETURNED", 0)

        if pairs_found == 0:
            st.warning("没有找到符合当前参数的引物。可以尝试：延长输入序列、降低 GC/Tm 限制，或检查模式是否选择正确。")
            explain = results.get("PRIMER_PAIR_EXPLAIN", "")
            if explain:
                st.caption(f"Primer3 反馈：{explain}")
        else:
            st.success(f"成功筛选出 {pairs_found} 对候选引物。")
            st.caption(status_msg)

            for i in range(pairs_found):
                f_seq = results[f"PRIMER_LEFT_{i}_SEQUENCE"]
                r_seq = results[f"PRIMER_RIGHT_{i}_SEQUENCE"]
                # 优先使用 Primer3 返回的 GC_PERCENT；如果某些版本没有该键，再用手动计算兜底。
                f_gc = results.get(f"PRIMER_LEFT_{i}_GC_PERCENT", calculate_gc(f_seq))
                r_gc = results.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", calculate_gc(r_seq))
                f_tm = results[f"PRIMER_LEFT_{i}_TM"]
                r_tm = results[f"PRIMER_RIGHT_{i}_TM"]
                f_gc_count = f_seq.count("G") + f_seq.count("C")
                r_gc_count = r_seq.count("G") + r_seq.count("C")
                prod_size = results[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]

                with st.expander(f"方案 {i + 1} (产物: {prod_size} bp)", expanded=(i == 0)):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("**Forward Primer**")
                        st.code(f_seq, language="text")
                        st.text(f"Length: {len(f_seq)} nt | Tm: {f_tm:.2f}°C | GC: {f_gc:.2f}% ({f_gc_count}/{len(f_seq)})")
                    with col2:
                        st.markdown("**Reverse Primer**")
                        st.code(r_seq, language="text")
                        st.text(f"Length: {len(r_seq)} nt | Tm: {r_tm:.2f}°C | GC: {r_gc:.2f}% ({r_gc_count}/{len(r_seq)})")

    except Exception as e:
        st.error(f"运行出错：{e}")
        st.caption("建议先确认：输入序列是否足够长、是否选择了正确模式、服务器环境是否已安装 primer3-py。")
