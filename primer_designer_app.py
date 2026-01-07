import streamlit as st
import primer3

# --- 页面配置 ---
st.set_page_config(page_title="Chanz帮你一秒设计引物", page_icon="🧬")

st.title("Chanz觉得你的时间很宝贵")
st.markdown("Written by Chanz | 模式 1 & 2 专为 circRNA 优化 ｜ 常规引物设计请选择模式3！")

# --- 侧边栏：参数设置 ---
with st.sidebar:
    st.header("⚙️ 参数设置")
    mode = st.radio(
        "选择设计模式",
        (1, 2, 3),
        format_func=lambda x: {
            1: "1. circRNA Divergent (qPCR)",
            2: "2. circRNA 跨接头 (Sanger测序)",
            3: "3. 常规线性 mRNA (普通qPCR)"
        }[x]
    )
    st.info("模式说明：\n\n**模式1**: 模拟环化，产物 80-150bp\n\n**模式2**: 模拟环化，产物 250-600bp\n\n**模式3**: 常规线性设计")

# --- 主界面 ---
raw_input = st.text_area("在此粘贴基因序列 (5'->3')", height=150, placeholder="ATGC...")
clean_seq = raw_input.replace(" ", "").replace("\n", "").replace("\r", "").upper()

if st.button("🚀 开始设计引物"):
    if not clean_seq:
        st.error("请输入序列！")
    else:
        seq_len = len(clean_seq)
        # 初始化基础参数
        params = {
            'SEQUENCE_ID': 'Target',
            'PRIMER_NUM_RETURN': 5,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_MAX_POLY_X': 5,
        }
        
        final_sequence = clean_seq
        status_msg = ""

        # 模式逻辑
        if mode == 1:
            final_sequence = clean_seq + clean_seq
            join_point = seq_len
            params['SEQUENCE_TARGET'] = [join_point - 10, 20]
            params['PRIMER_PRODUCT_SIZE_RANGE'] = [[80, 150]]
            status_msg = f"模式1：已执行序列串联模拟环化"
            
        elif mode == 2:
            final_sequence = clean_seq + clean_seq
            join_point = seq_len
            params['SEQUENCE_TARGET'] = [join_point - 5, 10]
            params['PRIMER_PRODUCT_SIZE_RANGE'] = [[250, 600], [200, 800]]
            status_msg = f"模式2：Sanger测序模式 (长片段)"

        elif mode == 3:
            params['PRIMER_PRODUCT_SIZE_RANGE'] = [[80, 200]]
            status_msg = "模式3：常规线性引物设计"

        params['SEQUENCE_TEMPLATE'] = final_sequence

        try:
            results = primer3.bindings.design_primers(
                {'SEQUENCE_ID': 'Target', 'SEQUENCE_TEMPLATE': final_sequence},
                params
            )
            
            pairs_found = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            
            if pairs_found == 0:
                st.warning("未找到合适的引物，请尝试放宽条件或检查序列长度。")
            else:
                st.success(f"成功找到 {pairs_found} 对引物！")
                st.caption(status_msg)
                
                for i in range(pairs_found):
                    with st.expander(f"方案 {i+1} (产物: {results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']} bp)", expanded=(i==0)):
                        c1, c2 = st.columns(2)
                        with c1:
                            st.markdown("**Forward Primer**")
                            st.code(results[f'PRIMER_LEFT_{i}_SEQUENCE'])
                            st.text(f"Tm: {results[f'PRIMER_LEFT_{i}_TM']:.1f}°C | GC: {results[f'PRIMER_LEFT_{i}_GC_PERCENT']:.1f}%")
                        with c2:
                            st.markdown("**Reverse Primer**")
                            st.code(results[f'PRIMER_RIGHT_{i}_SEQUENCE'])
                            st.text(f"Tm: {results[f'PRIMER_RIGHT_{i}_TM']:.1f}°C | GC: {results[f'PRIMER_RIGHT_{i}_GC_PERCENT']:.1f}%")

        except Exception as e:
            st.error(f"运行出错: {e}")
