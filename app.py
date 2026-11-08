from dataclasses import asdict
from io import StringIO
import csv
import numpy as np
from matplotlib.figure import Figure
from shiny import App, reactive, render, ui
from model import Conditions, calculate, chromatogram

CSS = '''
body{background:#f5f3ed;color:#172d39;font-family:system-ui,sans-serif}
.container-fluid{max-width:1500px;margin:auto;padding:28px}
h1{font-size:clamp(30px,4vw,52px);font-weight:700;letter-spacing:-2px;margin:8px 0 14px}
h2,h3{font-weight:650} .eyebrow{font-size:12px;letter-spacing:2px;color:#007b78;font-weight:750}
.lead{max-width:800px;color:#51616a;font-size:18px}.card{border:1px solid #dfE3df;border-radius:14px;box-shadow:none}
.card-header{background:white;font-weight:650;padding:16px 20px}.sidebar{background:#fff}
.btn-primary{background:#007b78;border-color:#007b78}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:20px 0}
.metric{background:#172d39;color:white;border-radius:12px;padding:18px}.metric strong{display:block;font-size:29px}.metric span{font-size:12px;color:#c3d7dc}
.note{padding:16px;background:#e9f3ef;border-left:3px solid #00857d;border-radius:4px;margin:12px 0}
.equation{font:600 21px Georgia,serif;padding:16px;background:#f5f3ed;border-radius:8px;margin:12px 0}
.lesson{padding:24px;max-width:950px}.lesson p,.lesson li{line-height:1.7}.small{color:#596c76;font-size:13px}
@media(max-width:700px){.metrics{grid-template-columns:repeat(2,1fr)}.container-fluid{padding:14px}.metric strong{font-size:23px}}
'''

def number(id, label, value, low, high, step):
    return ui.input_numeric(id, label, value, min=low, max=high, step=step)

app_ui = ui.page_fluid(
    ui.tags.style(CSS),
    ui.div('CHROMATOGRAPHY / AN INTERACTIVE TEACHING LAB', class_='eyebrow'),
    ui.h1('How fast should a column run?'),
    ui.p('Explore the van Deemter equation. Watch band broadening change a separation, then find the trade-off between resolution and analysis time.', class_='lead'),
    ui.navset_tab(
        ui.nav_panel('Explore the separation',
            ui.layout_sidebar(
                ui.sidebar(
                    ui.h4('Set the flow'),
                    ui.input_slider('u','Linear velocity u (mm/s)',0.1,5,1,step=0.05),
                    ui.input_action_button('optimal','Use optimum velocity',class_='btn-primary'),
                    ui.h5('Column & analytes'),
                    ui.input_slider('L','Column length L (mm)',50,300,150,step=10),
                    ui.input_slider('alpha','Selectivity α = k₂ / k₁',1,1.3,1.12,step=0.005),
                    ui.input_slider('k','Retention factor k₁',0.5,10,3,step=0.1),
                    ui.accordion(ui.accordion_panel('Band-broadening coefficients',
                        number('A','A · multiple paths (µm)',8,0,30,1),
                        number('B','B · longitudinal diffusion (µm·mm/s)',6,0.5,20,0.5),
                        number('C','C · mass transfer (µm·s/mm)',6,1,20,0.5)),open=False),
                    ui.input_action_button('save','Save current comparison'),
                    ui.input_action_button('reset','Reset experiment'),
                    ui.download_button('download','Download current data'),
                    ui.p('Illustrative packed-column model; not an instrument method predictor.',class_='small'),width=310),
                ui.output_ui('metrics'),
                ui.output_ui('interpretation'),
                ui.layout_columns(
                    ui.card(ui.card_header('01 / Where does broadening come from?'),ui.output_plot('deemter',height='330px')),
                    ui.card(ui.card_header('02 / What happens to the separation?'),
                            ui.input_switch('normalized','Normalize time by t₀ (compare peak overlap)',True),
                            ui.output_plot('peaks',height='290px')),
                    col_widths=(6,6)),
                ui.card(ui.card_header('03 / Resolution versus speed'),ui.output_plot('tradeoff',height='260px')),
                ui.output_ui('comparison')
            )),
        ui.nav_panel('Understand the model',ui.div(
            ui.h2('Three mechanisms. One efficiency curve.'),
            ui.div('H = A + B/u + Cu',class_='equation'),
            ui.p('H is the height equivalent to a theoretical plate. Smaller H means less band spreading per unit column length. A theoretical plate is a model of efficiency, not a physical disc inside the column.'),
            ui.h4('A · Multiple paths'),ui.p('Molecules take different routes through a packed bed. In this simplified model, this contribution is independent of velocity. Better packing can reduce it.'),
            ui.h4('B/u · Longitudinal diffusion'),ui.p('Molecules diffuse along the column axis. Slow flow gives diffusion more time to spread a band, so this term decreases as velocity increases.'),
            ui.h4('Cu · Resistance to mass transfer'),ui.p('Solutes need time to equilibrate between flowing mobile phase and less accessible regions or stationary phase. At high velocity, incomplete equilibration increases spreading.'),
            ui.div('uₒₚₜ = √(B/C)     •     Hₘᵢₙ = A + 2√(BC)     •     N = L/H',class_='equation'),
            ui.p('Use the same length units for L and H. The app converts millimetres to micrometres before calculating plate count N. At the optimum, B/u = Cu. Changing A changes the minimum height but does not shift the optimum velocity.'),
            ui.h2('From efficiency to resolution'),
            ui.div('t₀ = L/u     •     tR,i = t₀(1 + kᵢ)     •     σᵢ = tR,i/√N',class_='equation'),
            ui.div('Rs = 2(tR,2 − tR,1)/(w₁ + w₂), with wᵢ = 4σᵢ',class_='equation'),
            ui.p('The displayed Rs uses these Gaussian tangent-baseline widths, which keeps the number consistent with the simulated peaks. Gaussian tails never reach an exact zero; Rs ≈ 1.5 is a conventional practical baseline-separation target for comparable symmetric peaks.'),
            ui.div('Rs ≈ (√N/4) × [(α − 1)/α] × [k₂/(1 + k₂)]',class_='equation'),
            ui.p('This familiar teaching approximation uses the later peak’s width for both peaks. The app reports its value alongside the exact model calculation. With retention and selectivity fixed, Rs scales as √N: doubling column length improves resolution by √2, while doubling retention time at fixed velocity.'),
            ui.h4('What is held constant?'),
            ui.p('Isocratic conditions, Gaussian peaks of equal area, equal N for both analytes, fixed A/B/C, and no extra-column dispersion, overload, noise, pressure limit or temperature change. Retention and selectivity are independent controls here; real solvent and temperature changes can also change viscosity, diffusion and the coefficients.'),
            ui.p('u is linear mobile-phase velocity, not volumetric flow in mL/min. Converting requires column cross-sectional area and an appropriate mobile-phase porosity. Faster peaks can be narrower in seconds while still having worse resolution; use normalized time to see their relative overlap.'),
            ui.h4('Further reading'),
            ui.tags.a('Agilent · The LC Handbook',href='https://www.agilent.com/cs/library/primers/public/LC-Handbook-Complete-2.pdf',target='_blank'),
            class_='lesson')),
        ui.nav_panel('Guided experiments',ui.div(
            ui.h2('Predict. Change one variable. Explain.'),
            ui.p('Use Reset experiment before each activity. Save a comparison before changing a control.'),
            ui.h4('1 / Find the best flow'),ui.p('Compare u = 0.2, 1.0 and 4.0 mm/s. Which term dominates each extreme? Find the minimum H, then compare the resolution and retention time. Is the fastest run the best resolved?'),
            ui.h4('2 / Can a longer column rescue a separation?'),ui.p('Save the default 150 mm column, then double L to 300 mm. Predict the changes in N, Rs and retention time before reading them. Explain why resolution does not double.'),
            ui.h4('3 / When efficiency cannot help'),ui.p('Set α = 1. Move flow to the optimum and increase column length. Why are there still no separated peaks? Restore α = 1.12, then try 1.20 and compare this gain with doubling column length.'),
            ui.h4('4 / Move the optimum'),ui.p('Set B = 12 while keeping C = 6. Predict uₒₚₜ. Then increase A alone. Which edit moves the optimum, and which only lifts the curve?'),
            ui.h4('5 / Choose a practical method'),ui.p('At the default settings, find the highest velocity that still gives Rs ≥ 1.5. Compare its retention time with the minimum-H setting. Explain why best efficiency and shortest acceptable run are different objectives.'),
            ui.accordion(ui.accordion_panel('Instructor notes / expected observations',
                ui.p('Default uₒₚₜ = 1.00 mm/s; H = 20 µm; N = 7,500; exact Rs ≈ 1.865; tR,2 = 10.90 min. Doubling L multiplies N and time by 2 and Rs by √2. At α = 1, Rs = 0 for every efficiency. For B = 12 and C = 6, uₒₚₜ = √2 ≈ 1.414 mm/s. The default continuous-model upper velocity satisfying Rs ≥ 1.5 is about 3.53 mm/s; on the slider, 3.50 mm/s meets the target. This ignores pressure and other practical limits.')),open=False),
            class_='lesson')),
    ), title='Van Deemter Lab'
)

def server(input, output, session):
    reference = reactive.Value(None)

    @reactive.calc
    def settings():
        vals = {key:getattr(input,key)() for key in asdict(Conditions())}
        from shiny import req
        req(all(v is not None for v in vals.values()))
        req(0.1 <= vals['u'] <= 5 and 0 <= vals['A'] <= 30 and 0.5 <= vals['B'] <= 20 and 1 <= vals['C'] <= 20)
        return Conditions(**vals)

    @reactive.effect
    @reactive.event(input.optimal)
    def _optimal():
        p=settings()
        ui.update_slider('u',value=round((p.B/p.C)**0.5,4))

    @reactive.effect
    @reactive.event(input.save)
    def _save(): reference.set(settings())

    @reactive.effect
    @reactive.event(input.reset)
    def _reset():
        for key,value in asdict(Conditions()).items():
            (ui.update_numeric if key in ['A','B','C'] else ui.update_slider)(key,value=value)
        reference.set(None)

    @render.ui
    def metrics():
        m=calculate(settings())
        pairs=[('Plate height H',f"{m['H']:.1f} µm"),('Theoretical plates N',f"{m['N']:,.0f}"),('Resolution Rs',f"{m['Rs']:.2f}"),('Later peak tR,2',f"{m['t2']/60:.2f} min")]
        return ui.div(*[ui.div(ui.span(k),ui.tags.strong(v),class_='metric') for k,v in pairs],class_='metrics')

    @render.ui
    def interpretation():
        p=settings();m=calculate(p)
        state='Meets the Rs ≥ 1.5 target.' if m['Rs']>=1.5 else 'Below the Rs = 1.5 separation target.'
        terms={'multiple paths (A)':p.A,'longitudinal diffusion (B/u)':p.B/p.u,'mass transfer (Cu)':p.C*p.u}
        return ui.div(f"{state} Largest broadening contribution: {max(terms,key=terms.get)}. Optimum velocity: {m['uopt']:.2f} mm/s. Teaching approximation: Rs ≈ {m['Rs_approx']:.2f}.",class_='note')

    def axes():
        fig=Figure(figsize=(7,3.3),layout='constrained',facecolor='white');ax=fig.subplots()
        for side in ['top','right']: ax.spines[side].set_visible(False)
        ax.spines['left'].set_color('#d3dcdf');ax.spines['bottom'].set_color('#d3dcdf')
        ax.tick_params(colors='#51616a',labelsize=9);ax.grid(axis='y',alpha=.14)
        return fig,ax

    @render.plot(alt='Van Deemter curve with individual A, B/u and Cu contributions, current velocity and optimum.')
    def deemter():
        p=settings();m=calculate(p);fig,ax=axes();u=np.linspace(.1,5,400)
        for y,label,color,style in [(np.full_like(u,p.A),'A','#809198',':'),(p.B/u,'B/u','#4b78b8','--'),(p.C*u,'Cu','#c57a29','--'),(p.A+p.B/u+p.C*u,'Total H','#007b78','-')]:
            ax.plot(u,y,label=label,color=color,ls=style,lw=2)
        ax.axvline(m['uopt'],color='#8d9b9f',ls=':',label='Optimum')
        ax.scatter([p.u],[m['H']],color='#172d39',s=65,zorder=5,label='Current')
        ax.set(xlabel='Linear velocity u (mm/s)',ylabel='Plate height H (µm)',xlim=(.1,5),ylim=(0,None))
        ax.legend(ncol=3,fontsize=8,frameon=False);return fig

    @render.plot(alt='Two Gaussian chromatographic peaks and their summed detector signal; saved reference is dashed.')
    def peaks():
        p=settings();fig,ax=axes();norm=input.normalized()
        x,y1,y2=chromatogram(p,norm)
        ax.fill_between(x,y1,color='#007b78',alpha=.18);ax.plot(x,y1,color='#007b78',label='Solute 1')
        ax.fill_between(x,y2,color='#c57a29',alpha=.18);ax.plot(x,y2,color='#c57a29',label='Solute 2')
        ax.plot(x,y1+y2,color='#172d39',lw=1.5,label='Total signal')
        if reference() is not None:
            xr,a,b=chromatogram(reference(),norm);ax.plot(xr,a+b,color='#8a6ca7',ls='--',label='Saved total')
        ax.set(xlabel='Normalized time t/t₀' if norm else 'Time (min)',ylabel='Area-normalized signal',ylim=(0,None))
        ax.legend(fontsize=8,frameon=False);return fig

    @render.plot(alt='Resolution versus velocity with a 1.5 target and later-peak retention time on a second panel.')
    def tradeoff():
        p=settings();fig=Figure(figsize=(12,2.7),layout='constrained');axs=fig.subplots(1,2)
        u=np.linspace(.1,5,400);ms=[calculate(Conditions(**{**asdict(p),'u':v})) for v in u]
        for ax,key,label,color in zip(axs,['Rs','t2'],['Resolution Rs','Later-peak time (min)'],['#007b78','#c57a29']):
            scale=60 if key=='t2' else 1
            ax.plot(u,[m[key]/scale for m in ms],color=color,lw=2)
            ax.scatter([p.u],[calculate(p)[key]/scale],color='#172d39',zorder=5)
            ax.set(xlabel='Linear velocity u (mm/s)',ylabel=label,xlim=(.1,5));ax.grid(alpha=.15)
            ax.spines[['top','right']].set_visible(False)
        axs[0].axhline(1.5,color='#c57a29',ls='--',label='Rs = 1.5');axs[0].legend(frameon=False,fontsize=8)
        return fig

    @render.ui
    def comparison():
        old=reference()
        if old is None:return ui.p('Save a comparison to measure what changed.',class_='small')
        a,b=calculate(old),calculate(settings())
        rows=[ui.tags.tr(ui.tags.th('Measurement'),ui.tags.th('Saved'),ui.tags.th('Current'))]
        for label,x,y in [('u (mm/s)',old.u,settings().u),('H (µm)',a['H'],b['H']),('N',a['N'],b['N']),('Rs',a['Rs'],b['Rs']),('tR,2 (min)',a['t2']/60,b['t2']/60)]:
            rows.append(ui.tags.tr(ui.tags.td(label),ui.tags.td(f'{x:.3f}'),ui.tags.td(f'{y:.3f}')))
        return ui.card(ui.card_header('Saved comparison'),ui.tags.table(*rows,class_='table'))

    @render.download(filename='van-deemter-experiment.csv')
    def download():
        p=settings();m=calculate(p);s=StringIO();w=csv.writer(s)
        w.writerow(['parameter','value']);w.writerows(asdict(p).items());w.writerows(m.items())
        w.writerow([]);w.writerow(['time_min','solute_1_area_normalized','solute_2_area_normalized','total'])
        x,a,b=chromatogram(p)
        w.writerows(zip(x,a,b,a+b));yield s.getvalue()

app=App(app_ui,server)
