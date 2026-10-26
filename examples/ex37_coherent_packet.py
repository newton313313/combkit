"""
範例 37：相干態波包在諧振子裡來回不散（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 8 章 §8.2–8.3、
(8.2.15)–(8.2.28)、(8.3.1)–(8.3.9)、Fig. 8.2.1–8.3.2）——先用程式算：|α|² = 5（x0 = √10、p0 = 0，ħ = μ = ω = 1）的 Poisson 振幅
c_n = e^{−|α|²/2}αⁿ/√n!（Q|0⟩ 的 90 維截斷矩陣 ＝ 閉式），⟨n⟩ = 5、Δn = √5；波包 χ(x, t) 的閉式 (8.2.28) ＝ 本徵函數的和 ＝ FFT 分裂算子的數值演化；
十張等間隔快照乘上 e^{i(5+½)t_p} 相加就得到 ⟨x|5⟩，誤差 c15/c5 = 0.0299。畫：相空間的相量順時針轉一圈（PhasePhasor）、
波包的 Re χ 與 |χ| 跟著來回、右邊 Poisson 長條；最後十張快照的和與 ⟨x|5⟩ 疊在一起。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex37_coherent_packet.py CoherentPacketExample
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

import numpy as np
from manimlib import *

from combkit import sym_dyn as DY
from combkit.draw import sym_draw as SD
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene

XS = np.linspace(-6.0, 6.0, 110)


def check():
    """c_n（閉式 ＝ Q|0⟩）、⟨n⟩ 與 Δn、波包（閉式 ＝ 本徵函數和 ＝ 分裂算子）、快照合成的誤差。回傳 (x0, c, 快照和, ⟨x|5⟩)。"""
    x0 = DY.coh_x0()
    al = DY.alpha_of(x0, 0.0)
    c = DY.coherent_cn(al, 15)
    assert np.allclose(DY.coherent_numeric(x0, 0.0)[:16], c, atol=1e-9)
    st = DY.coherent_stats(al)
    assert abs(st["mean"] - 5) < 1e-10 and abs(st["dn"] - math.sqrt(5)) < 1e-10
    assert f"{abs(c[4]):.4f}" == f"{abs(c[5]):.4f}" == "0.4189"
    for t in (0.0, 1.3, math.pi):
        assert np.allclose(DY.packet_closed(XS, t, x0), DY.packet_sum(XS, t, x0), atol=1e-9)
    xg, ps = DY.packet_split_operator(2.0, x0)
    assert np.allclose(ps, DY.packet_closed(xg, 2.0, x0), atol=5e-5)
    se = DY.snapshot_error(x0, 5, 10)
    assert f"{se['ratio']:.4f}" == "0.0299" and abs(se["l2"] - se["pred"]) < 1e-6
    return x0, c, DY.snapshot_sum(XS, x0, 5, 10).real, (c[5] * DY.ho_eigfun(5, XS)).real


def _packet_updater(tr, x0, ux, uy, base):
    def upd(m):
        psi = DY.packet_closed(XS, tr.get_value(), x0)
        m[0].set_points_as_corners([base + np.array([x * ux, v * uy, 0]) for x, v in zip(XS, psi.real)])
        m[1].set_points_as_corners([base + np.array([x * ux, v * uy, 0]) for x, v in zip(XS, abs(psi))])
    return upd


def _phasor_updater(tr, amp):
    def upd(m):
        m.set_alpha(amp * np.exp(-1j * tr.get_value()))
    return upd


class CoherentPacketExample(ChuScene):
    def construct(self):
        # 1. 數學
        x0, c, ssum, exact = check()
        amp = math.sqrt(5)

        # 2. 畫
        self.show_header("相干態：來回不散的波包")
        ph = SD.PhasePhasor(radius=1.1, center=LEFT * 4.6 + UP * 0.8, amp=amp)
        base = np.array([0.2, 0.3, 0.0])
        ux, uy = 0.42, 1.6
        psi0 = DY.packet_closed(XS, 0.0, x0)
        wave = VGroup(VMobject().set_points_as_corners([base + np.array([x * ux, v * uy, 0]) for x, v in zip(XS, psi0.real)])
                      .set_stroke(SD.PACKET_COLOR, 2.5),
                      VMobject().set_points_as_corners([base + np.array([x * ux, v * uy, 0]) for x, v in zip(XS, abs(psi0))])
                      .set_stroke(SD.PACKET_ENV_COLOR, 1.5))
        axis = Line(base + LEFT * 6 * ux, base + RIGHT * 6 * ux).set_stroke(GREY_D, 1.5)
        bars = SD.ProbBars([abs(v) for v in c[:11]], labels=[str(n) for n in range(11)], width=0.2, gap=0.06, height=1.8,
                           font_size=16, label_every=5).move_to(RIGHT * 4.9 + DOWN * 1.6)
        for m in (ph, wave, bars):
            assert_in_frame(m)
        self.play(FadeIn(ph), ShowCreation(axis), ShowCreation(wave), run_time=1.2)
        self.say("$|\\alpha|^2=5$", "：基態被拉到 ", "$x_0=\\sqrt{10}$", "，放手")
        tr = ValueTracker(0.0)
        wave.add_updater(_packet_updater(tr, x0, ux, uy, base))
        ph.add_updater(_phasor_updater(tr, amp))
        self.play(tr.animate.set_value(2 * math.pi), run_time=6.0, rate_func=linear)
        wave.clear_updaters()
        ph.clear_updaters()
        self.say("相量 ", "$\\alpha_t=\\alpha e^{-it}$", " 順時針轉，包絡形狀不變")
        self.play(FadeIn(bars), run_time=1.0)
        self.say("Poisson：", "$c_4=c_5=0.4189$", "，", "$\\Delta n=\\sqrt5$")
        crv_s = VMobject().set_points_as_corners([base + DOWN * 2.3 + np.array([x * ux, v * 2.4, 0]) for x, v in zip(XS, ssum)]).set_stroke(SD.PACKET_COLOR, 3)
        crv_e = DashedVMobject(VMobject().set_points_as_corners([base + DOWN * 2.3 + np.array([x * ux, v * 2.4, 0]) for x, v in zip(XS, exact)]),
                               num_dashes=50).set_stroke(YELLOW, 2)
        assert_in_frame(VGroup(crv_s, crv_e))
        self.play(ShowCreation(crv_s), ShowCreation(crv_e), run_time=1.5)
        self.say("十張快照乘相位相加 ＝ ", "$\\langle x|5\\rangle$", "，誤差 ", "$c_{15}/c_5=0.0299$")
        self.wait(1.0)
