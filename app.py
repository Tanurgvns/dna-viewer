"""DNA Viewer: a small Streamlit app for exploring DNA sequences visually."""

import random

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import dna_utils as du
from helix_view import helix_html

BASE_COLOURS = {
    "A": "#3DDC84",  # green
    "T": "#4FC3F7",  # blue
    "G": "#FFCA28",  # yellow
    "C": "#B388FF",  # purple
    "N": "#9E9E9E",  # grey
}
MAX_DISPLAY = 6000  # bases shown at once in the coloured view
MAX_HELIX = 1_000_000  # bases sent to the helix view

st.set_page_config(page_title="DNA Viewer", page_icon="🧬", layout="wide")

st.markdown(
    """
    <style>
    .stApp {background: #000;}
    .block-container {padding-top: 3.5rem;}
    .seq-block {font-family: "JetBrains Mono", "Fira Code", Consolas, monospace;
                font-size: 14px; line-height: 1.7; overflow-x: auto;
                white-space: pre; padding: 0.75rem 1rem; border-radius: 6px;
                border: 1px solid rgba(128,128,128,0.25);}
    .seq-block .pos {color: #9AA0A6; user-select: none;}
    .seq-block .hit {outline: 2px solid currentColor; border-radius: 2px;}
    .legend span {display: inline-block; margin-right: 1rem; font-weight: 600;}
    </style>
    """,
    unsafe_allow_html=True,
)


# ---------- helpers ----------

@st.cache_data
def example_fasta() -> str:
    """A made-up demo sequence with a couple of ORFs planted in it."""
    rng = random.Random(42)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    sense_codons = [c for c, aa in du.CODON_TABLE.items() if aa != "*"]
    orf = lambda n: "ATG" + "".join(rng.choice(sense_codons) for _ in range(n)) + "TAA"
    seq = rand(180) + orf(120) + rand(250) + du.reverse_complement(orf(80)) + rand(200)
    seq = seq[:420] + "GAATTC" + seq[420:]  # an EcoRI site to search for
    lines = [seq[i:i + 70] for i in range(0, len(seq), 70)]
    return ">demo_sequence synthetic example for DNA Viewer\n" + "\n".join(lines)


def coloured_sequence_html(seq: str, offset: int, width: int,
                           highlight: set[int] | None = None,
                           show_complement: bool = False) -> str:
    highlight = highlight or set()
    rows = []
    pad = len(str(offset + len(seq)))
    for i in range(0, len(seq), width):
        chunk = seq[i:i + width]
        spans = []
        for j, base in enumerate(chunk):
            cls = ' class="hit"' if (offset + i + j) in highlight else ""
            spans.append(f'<span style="color:{BASE_COLOURS[base]}"{cls}>{base}</span>')
        rows.append(f'<span class="pos">{offset + i + 1:>{pad}}  </span>{"".join(spans)}')
        if show_complement:
            comp = "".join(
                f'<span style="color:{BASE_COLOURS[b]};opacity:0.55">{b}</span>'
                for b in du.complement(chunk)
            )
            rows.append(f'<span class="pos">{"":>{pad}}  </span>{comp}\n')
    return f'<div class="seq-block">{"<br>".join(rows)}</div>'


def legend_html() -> str:
    names = {"A": "Adenine", "T": "Thymine", "G": "Guanine", "C": "Cytosine", "N": "Unknown"}
    return '<div class="legend">' + "".join(
        f'<span style="color:{c}">{b} {names[b]}</span>' for b, c in BASE_COLOURS.items()
    ) + "</div>"


# ---------- sidebar: input ----------

with st.sidebar:
    st.title("🧬 DNA Viewer")
    source = st.radio("Sequence source", ["Example", "Upload FASTA", "Paste sequence"])

    text = ""
    file_name = "Pasted sequence"
    if source == "Example":
        text = example_fasta()
        file_name = "example.fasta"
    elif source == "Upload FASTA":
        uploaded = st.file_uploader("FASTA file", type=["fa", "fasta", "fna", "ffn", "txt"])
        if uploaded:
            text = uploaded.getvalue().decode("utf-8", errors="replace")
            file_name = uploaded.name
    else:
        text = st.text_area("Paste FASTA or raw sequence", height=200,
                            placeholder=">my_seq\nATGCGTACGTTAG...")

records = du.parse_fasta(text)

if not records:
    st.title("🧬 DNA Viewer")
    st.info("Choose the example, upload a FASTA file, or paste a sequence in the sidebar to start.")
    st.stop()

with st.sidebar:
    if len(records) > 1:
        idx = st.selectbox("Record", range(len(records)),
                           format_func=lambda i: f"{records[i].header[:50]} ({len(records[i].sequence):,} bp)")
    else:
        idx = 0
    st.caption(f"{len(records)} record(s) loaded")

record = records[idx]
seq = record.sequence

helix_seq = seq[:MAX_HELIX]
if len(seq) > MAX_HELIX:
    st.caption(f"The helix view shows the first {MAX_HELIX:,} bases.")
page = helix_html(helix_seq, file_name)
if hasattr(st, "iframe"):
    st.iframe(page, height=820)
else:  # older Streamlit versions
    import streamlit.components.v1 as components
    components.html(page, height=820, scrolling=True)

st.subheader("More analysis")
st.caption(record.header)

tab_overview, tab_seq, tab_gc, tab_orf, tab_trans, tab_motif = st.tabs(
    ["Overview", "Sequence", "GC content", "ORF finder", "Translation", "Motif search"]
)

# ---------- overview ----------

with tab_overview:
    counts = du.base_counts(seq)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Length", f"{len(seq):,} bp")
    c2.metric("GC content", f"{du.gc_content(seq):.1f}%")
    c3.metric("GC skew", f"{du.gc_skew(seq):+.3f}")
    c4.metric("Approx. dsDNA weight", f"{du.molecular_weight_dsdna(seq) / 1000:,.1f} kDa")

    if counts["N"]:
        st.warning(f"{counts['N']:,} unknown bases (N) in this sequence.")

    left, right = st.columns([3, 2])
    with left:
        df = pd.DataFrame({"Base": list(counts), "Count": list(counts.values())})
        df = df[df["Count"] > 0]
        fig = px.bar(df, x="Base", y="Count", color="Base",
                     color_discrete_map=BASE_COLOURS, text="Count")
        fig.update_layout(showlegend=False, height=340, margin=dict(t=30, b=10),
                          title="Base composition")
        st.plotly_chart(fig, width="stretch")
    with right:
        fig = px.pie(df, names="Base", values="Count", color="Base",
                     color_discrete_map=BASE_COLOURS, hole=0.55)
        fig.update_layout(height=340, margin=dict(t=30, b=10), title="Proportions")
        st.plotly_chart(fig, width="stretch")

# ---------- coloured sequence ----------

with tab_seq:
    st.markdown(legend_html(), unsafe_allow_html=True)
    o1, o2, o3 = st.columns(3)
    width = o1.select_slider("Bases per line", [30, 50, 60, 80, 100], value=60)
    strand = o2.radio("Strand", ["Forward", "Reverse complement"], horizontal=True)
    show_comp = o3.toggle("Show complement strand", value=False)

    view_seq = seq if strand == "Forward" else du.reverse_complement(seq)
    if len(view_seq) > MAX_DISPLAY:
        start = st.slider("Start position", 1, len(view_seq) - MAX_DISPLAY + 1, 1, step=width)
        st.caption(f"Showing {MAX_DISPLAY:,} of {len(view_seq):,} bases.")
    else:
        start = 1
    window = view_seq[start - 1:start - 1 + MAX_DISPLAY]
    st.markdown(coloured_sequence_html(window, start - 1, width, show_complement=show_comp),
                unsafe_allow_html=True)

    st.download_button("Download this strand as FASTA",
                       f">{record.header} ({strand.lower()})\n{view_seq}\n",
                       file_name="sequence.fasta")

# ---------- GC content ----------

with tab_gc:
    max_win = max(10, len(seq))
    g1, g2 = st.columns(2)
    win = g1.slider("Window size (bp)", 10, min(2000, max_win), min(100, max_win))
    step = g2.slider("Step (bp)", 1, max(1, win), max(1, win // 5))
    data = pd.DataFrame(du.sliding_window(seq, win, step),
                        columns=["Position", "GC %", "GC skew"])

    fig = go.Figure()
    fig.add_scatter(x=data["Position"], y=data["GC %"], mode="lines",
                    line=dict(color=BASE_COLOURS["G"], width=2), name="GC %")
    fig.add_hline(y=du.gc_content(seq), line_dash="dot",
                  annotation_text=f"Overall {du.gc_content(seq):.1f}%")
    fig.update_layout(height=320, title="GC content along the sequence",
                      xaxis_title="Position (bp)", yaxis_title="GC %",
                      margin=dict(t=40, b=10))
    st.plotly_chart(fig, width="stretch")

    fig = go.Figure()
    fig.add_scatter(x=data["Position"], y=data["GC skew"], mode="lines",
                    line=dict(color=BASE_COLOURS["C"], width=2), name="GC skew")
    fig.add_scatter(x=data["Position"], y=data["GC skew"].cumsum(), mode="lines",
                    line=dict(color=BASE_COLOURS["T"], width=1.5, dash="dash"),
                    name="Cumulative skew", yaxis="y2")
    fig.update_layout(height=320, title="GC skew (G−C)/(G+C)",
                      xaxis_title="Position (bp)", yaxis_title="Skew",
                      yaxis2=dict(overlaying="y", side="right", title="Cumulative"),
                      margin=dict(t=40, b=10), legend=dict(orientation="h", y=-0.25))
    st.plotly_chart(fig, width="stretch")
    st.caption("In bacterial genomes, the minimum and maximum of cumulative GC skew "
               "often mark the origin and terminus of replication.")

# ---------- ORFs ----------

with tab_orf:
    min_aa = st.slider("Minimum ORF length (amino acids)", 10, 300, 30)
    orfs = du.find_orfs(seq, min_aa)
    st.write(f"**{len(orfs)}** ORFs found (ATG start to stop codon, both strands).")

    if orfs:
        fig = go.Figure()
        fig.add_shape(type="line", x0=1, x1=len(seq), y0=0, y1=0,
                      line=dict(color="#9AA0A6", width=2))
        for o in orfs:
            y = o.frame if o.strand == "+" else -o.frame
            colour = BASE_COLOURS["A"] if o.strand == "+" else BASE_COLOURS["T"]
            fig.add_shape(type="rect", x0=o.start, x1=o.end, y0=y - 0.35, y1=y + 0.35,
                          fillcolor=colour, line=dict(width=0), opacity=0.85)
            fig.add_scatter(x=[(o.start + o.end) / 2], y=[y], mode="markers",
                            marker=dict(size=1, color=colour), showlegend=False,
                            hovertext=f"{o.strand}{o.frame}: {o.start}–{o.end} ({o.length_aa} aa)",
                            hoverinfo="text")
        fig.update_layout(
            height=320, title="ORF map",
            xaxis=dict(title="Position (bp)", range=[1, len(seq)]),
            yaxis=dict(tickvals=[-3, -2, -1, 0, 1, 2, 3],
                       ticktext=["−3", "−2", "−1", "", "+1", "+2", "+3"],
                       title="Strand / frame", range=[-3.8, 3.8]),
            margin=dict(t=40, b=10),
        )
        st.plotly_chart(fig, width="stretch")

        table = pd.DataFrame([{
            "Strand": o.strand, "Frame": o.frame, "Start": o.start, "End": o.end,
            "Length (nt)": o.length_nt, "Length (aa)": o.length_aa, "Protein": o.protein,
        } for o in orfs])
        st.dataframe(table, width="stretch", hide_index=True)
        fasta_out = "\n".join(f">ORF_{i + 1} {o.strand}{o.frame} {o.start}-{o.end}\n{o.protein}"
                              for i, o in enumerate(orfs))
        st.download_button("Download proteins as FASTA", fasta_out, file_name="orfs.faa")
    else:
        st.info("No ORFs at this length. Try lowering the minimum.")

# ---------- translation ----------

with tab_trans:
    st.caption("All six reading frames. * marks a stop codon, M a methionine (possible start).")
    frames = du.six_frame_translation(seq)
    for name, protein in frames.items():
        shown = protein[:2000]
        html = "".join(
            f'<span style="color:{BASE_COLOURS["T"]};font-weight:700">*</span>' if aa == "*"
            else f'<span style="color:{BASE_COLOURS["A"]};font-weight:700">M</span>' if aa == "M"
            else aa
            for aa in shown
        )
        st.markdown(f"**Frame {name}** ({len(protein):,} aa)")
        st.markdown(f'<div class="seq-block" style="white-space:pre-wrap;word-break:break-all">{html}</div>',
                    unsafe_allow_html=True)

# ---------- motif search ----------

with tab_motif:
    m1, m2 = st.columns([3, 1])
    motif = m1.text_input("Motif (IUPAC codes allowed, e.g. GAATTC or TATAWAW)", "GAATTC")
    both = m2.toggle("Search both strands", value=True)

    if motif:
        try:
            hits = du.find_motif(seq, motif, both)
        except ValueError as err:
            st.error(str(err))
            hits = None
        if hits is not None:
            st.write(f"**{len(hits)}** match(es) for `{motif.upper()}`.")
            if hits:
                st.dataframe(pd.DataFrame(hits), width="stretch", hide_index=True)
                highlight = {p - 1 for h in hits for p in range(h["start"], h["end"] + 1)}
                first = min(h["start"] for h in hits)
                ctx_start = max(0, first - 1 - 120)
                ctx = seq[ctx_start:ctx_start + 600]
                st.markdown("Context around the first match (matches outlined):")
                st.markdown(coloured_sequence_html(ctx, ctx_start, 60, highlight),
                            unsafe_allow_html=True)

st.divider()
st.caption("DNA Viewer is a small coding experiment, built with Streamlit.")
