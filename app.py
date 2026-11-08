from dataclasses import asdict
from io import StringIO
import csv
import numpy as np
from matplotlib.figure import Figure
from shiny import App, reactive, render, ui
from model import Conditions, calculate, chromatogram

CSS = '''
body{background:#f5f3ed;color:#202124;font-family:system-ui,sans-serif}
.container-fluid{max-width:1500px;margin:auto;padding:28px}
h1{font-size:clamp(30px,4vw,52px);font-weight:700;letter-spacing:-2px;margin:8px 0 14px}
h1,h2,h3{font-family:Georgia,serif} h2,h3{font-weight:650} .eyebrow{font-size:12px;letter-spacing:2px;color:#b5122a;font-weight:750}
.lead{max-width:800px;color:#51616a;font-size:18px}.card{border:1px solid #dfE3df;border-radius:14px;box-shadow:none}
.card-header{background:white;font-weight:650;padding:16px 20px}.sidebar{background:#fff}
.btn-primary{background:#b5122a;border-color:#b5122a}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:20px 0}
.metric{background:#202124;color:white;border-radius:12px;padding:18px}.metric strong{display:block;font-size:29px}.metric span{font-size:12px;color:#eee4df}
.note{padding:16px;background:#faf0ee;border-left:3px solid #b5122a;border-radius:4px;margin:12px 0}
.equation{font:600 21px Georgia,serif;padding:16px;background:#f5f3ed;border-radius:8px;margin:12px 0}
.controls{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}.control{padding:16px;border:1px solid #e8e0dc;border-radius:8px;background:#fff}.control p{font-size:13px;color:#626262;margin:6px 0}.irs--shiny .irs-bar,.irs--shiny .irs-single{background:#b5122a;border-color:#b5122a}.workbench{display:grid;grid-template-columns:minmax(300px,4fr) minmax(0,6fr);gap:20px;align-items:start}.toolbar{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0}.nav-link{color:#a10e25}.metric small{display:block;color:#eee4df;margin-top:6px}.lesson{padding:24px;max-width:950px}.lesson p,.lesson li{line-height:1.7}.small{color:#596c76;font-size:13px}
@media(max-width:900px){.workbench{grid-template-columns:1fr}} @media(max-width:700px){.controls{grid-template-columns:1fr}.metrics{grid-template-columns:repeat(2,1fr)}.container-fluid{padding:14px}.metric strong{font-size:23px}}
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
            ui.div(
                ui.div(
                    ui.input_select('experiment','Choose an investigation',
                        {'flow':'1 · Why not run as fast as possible?',
                         'length':'2 · Is a longer column always better?',
                         'selectivity':'3 · Can efficiency separate identical retention?',
                         'free':'4 · Design your own separation'}),
                    ui.input_action_button('load','Start investigation',class_='btn-primary'),
                    class_='toolbar'),
                ui.output_ui('prompt'),
                ui.div(
                    ui.card(ui.card_header('CHANGE ONE VARIABLE · Drag a slider'),
                        ui.div(
                            ui.div(ui.input_slider('u','Mobile-phase velocity · u (mm/s)',.1,5,1,step=.05),
                                ui.p('Slow ← → Fast. Predict: will faster flow sharpen or broaden the bands?'),class_='control'),
                            ui.div(ui.input_slider('L','Column length · L (mm)',50,300,150,step=10),
                                ui.p('Short ← → Long. More column gives more plates, but takes more time.'),class_='control'),
                            ui.div(ui.input_slider('alpha','Selectivity · α = k₂/k₁',1,1.3,1.12,step=.005),
                                ui.p('Similar ← → Different retention. In the lab, stationary phase or solvent composition can change selectivity.'),class_='control'),
                            ui.div(ui.input_slider('k','Retention factor · k₁',.5,10,3,step=.1),
                                ui.p('Weak ← → Strong retention. In the lab, mobile-phase strength is one way to change retention.'),class_='control'),class_='controls'),
                        ui.div(ui.input_action_button('optimal','Set best-efficiency flow'),
                               ui.input_action_button('save','Use current as reference'),
                               ui.input_action_button('reset','Reset all'),class_='toolbar'),
                        ui.p('Selectivity and retention are independent model controls. Real solvent changes may affect both, as well as band broadening.',class_='small')),
                    ui.card(ui.card_header('OBSERVE · Do the peaks overlap?'),
                        ui.input_switch('normalized','Compare overlap on normalized time (t/t₀)',True),
                        ui.output_plot('peaks',height='340px'),
                        ui.p('Solid curves: current experiment. Dashed curve: reference. Each solute has the same total peak area. Switch to minutes to see the analysis time.',class_='small')),
                    class_='workbench'),
                ui.output_ui('metrics'),ui.output_ui('interpretation'),
                ui.layout_columns(
                    ui.card(ui.card_header('EXPLAIN · Why does the band spread?'),ui.output_plot('deemter',height='330px')),
                    ui.card(ui.card_header('CONNECT · From broadening to separation'),
                        ui.output_ui('worked'),
                        ui.p('Smaller H → more plates N → less relative peak width → greater resolution, if retention and selectivity stay fixed.'),
                        ui.accordion(ui.accordion_panel('Explore the three broadening mechanisms',
                            ui.input_slider('A','Multiple paths · A (µm)',0,30,8,step=1),
                            ui.p('Different routes through the packing. Changes H without shifting the optimum.'),
                            ui.input_slider('B','Longitudinal diffusion · B (µm·mm/s)',.5,20,6,step=.5),
                            ui.p('More time at slow flow allows more diffusion along the column.'),
                            ui.input_slider('C','Mass-transfer resistance · C (µm·s/mm)',1,20,6,step=.5),
                            ui.p('At fast flow, solutes have less time to equilibrate. These are model coefficients, not instrument settings.')),open=False)),
                    col_widths=(6,6)),
                ui.card(ui.card_header('TRADE-OFF · Separation quality and time'),ui.output_plot('tradeoff',height='260px')),
                ui.output_ui('comparison'),ui.download_button('download','Download experiment data')
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
    reference = reactive.Value(Conditions())

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
            ui.update_slider(key,value=value)
        reference.set(Conditions())

    @reactive.effect
    @reactive.event(input.load)
    def _load():
        choices={'flow':Conditions(u=.2),'length':Conditions(L=100),
                 'selectivity':Conditions(alpha=1),'free':Conditions()}
        p=choices[input.experiment()]
        for key,value in asdict(p).items():ui.update_slider(key,value=value)
        reference.set(p)

    @render.ui
    def prompt():
        prompts={
            'flow':('Predict before you drag','Start at slow flow, then increase u. Does resolution keep improving? Find where the diffusion and mass-transfer terms balance. Hold the other sliders fixed.'),
            'length':('Make a quantitative prediction','Start with a 100 mm column. Double its length. Predict the factors by which plate count, resolution and retention time change; compare with the reference.'),
            'selectivity':('Try to separate the inseparable','Start with α = 1. Try changing velocity and length first. Then change α. Explain why reducing peak width cannot separate peaks with identical retention times.'),
            'free':('Choose a defensible compromise','Aim for Rs ≥ 1.5 with a short analysis time. Save a reference, change one parameter, and explain both the benefit and cost. Pressure limits are not modeled.')}
        title,body=prompts[input.experiment()]
        return ui.div(ui.tags.strong(title),ui.p(body),class_='note')

    @render.ui
    def worked():
        p=settings();m=calculate(p)
        return ui.div(
            ui.div(f"H = {p.A:g} + {p.B:g}/{p.u:g} + {p.C:g} × {p.u:g} = {m['H']:.1f} µm",class_='equation'),
            ui.div(f"N = {p.L*1000:,.0f} µm / {m['H']:.1f} µm = {m['N']:,.0f}",class_='equation'),
            ui.p(f"Peak spacing: {(m['t2']-m['t1'])/60:.3f} min. Average tangent-baseline width: {(2*(m['s1']+m['s2']))/60:.3f} min. Their ratio gives Rs = {m['Rs']:.2f}."))

    @render.ui
    def metrics():
        m=calculate(settings());r=calculate(reference())
        pairs=[('Plate height H','H','µm',1),('Theoretical plates N','N','',1),('Resolution Rs','Rs','',1),('Later peak tR,2','t2','min',60)]
        cards=[]
        for label,key,unit,scale in pairs:
            val=m[key]/scale
            delta=f"{(m[key]/r[key]-1)*100:+.0f}% vs reference" if r[key] else f"Reference: {r[key]:.2f}"
            cards.append(ui.div(ui.span(label),ui.tags.strong(f"{val:,.0f} {unit}" if key=='N' else f"{val:.2f} {unit}"),ui.tags.small(delta),class_='metric'))
        return ui.div(*cards,class_='metrics')

    @render.ui
    def interpretation():
        p=settings();m=calculate(p)
        state='Meets the Rs ≥ 1.5 target.' if m['Rs']>=1.5 else 'Below the Rs = 1.5 separation target.'
        terms={'multiple paths (A)':p.A,'longitudinal diffusion (B/u)':p.B/p.u,'mass transfer (Cu)':p.C*p.u}
        regime = ('The solutes have identical retention: efficiency alone cannot separate them.' if p.alpha == 1 else
                  'Below the optimum: increasing velocity reduces H and improves resolution.' if p.u < m['uopt']*.98 else
                  'Above the optimum: increasing velocity raises H and reduces resolution, even though the run is faster.' if p.u > m['uopt']*1.02 else
                  'Near the minimum plate height: longitudinal diffusion and mass transfer are balanced.')
        return ui.div(ui.tags.strong(regime),ui.p(f"{state} Largest broadening contribution: {max(terms,key=terms.get)}. Optimum velocity: {m['uopt']:.2f} mm/s. Teaching approximation: Rs ≈ {m['Rs_approx']:.2f}."),class_='note')

    def axes():
        fig=Figure(figsize=(7,3.3),layout='constrained',facecolor='white');ax=fig.subplots()
        for side in ['top','right']: ax.spines[side].set_visible(False)
        ax.spines['left'].set_color('#d3dcdf');ax.spines['bottom'].set_color('#d3dcdf')
        ax.tick_params(colors='#51616a',labelsize=9);ax.grid(axis='y',alpha=.14)
        return fig,ax

    @render.plot(alt='Van Deemter curve with individual A, B/u and Cu contributions, current velocity and optimum.')
    def deemter():
        p=settings();m=calculate(p);fig,ax=axes();u=np.linspace(.1,5,400)
        for y,label,color,style in [(np.full_like(u,p.A),'A','#809198',':'),(p.B/u,'B/u','#4b78b8','--'),(p.C*u,'Cu','#c57a29','--'),(p.A+p.B/u+p.C*u,'Total H','#b5122a','-')]:
            ax.plot(u,y,label=label,color=color,ls=style,lw=2)
        ax.axvline(m['uopt'],color='#8d9b9f',ls=':',label='Optimum')
        old=reference();ax.scatter([old.u],[calculate(old)['H']],marker='x',color='#8b8581',s=60,label='Reference',zorder=4)
        ax.scatter([p.u],[m['H']],color='#202124',s=65,zorder=5,label='Current')
        ax.set(xlabel='Linear velocity u (mm/s)',ylabel='Plate height H (µm)',xlim=(.1,5),ylim=(0,None))
        ax.legend(ncol=3,fontsize=8,frameon=False);return fig

    @render.plot(alt='Two Gaussian chromatographic peaks and their summed detector signal; saved reference is dashed.')
    def peaks():
        p=settings();fig,ax=axes();norm=input.normalized()
        x,y1,y2=chromatogram(p,norm)
        ax.fill_between(x,y1,color='#b5122a',alpha=.18);ax.plot(x,y1,color='#b5122a',label='Solute 1')
        ax.fill_between(x,y2,color='#c57a29',alpha=.18);ax.plot(x,y2,color='#c57a29',label='Solute 2')
        ax.plot(x,y1+y2,color='#202124',lw=1.5,label='Total signal')
        if reference() is not None:
            xr,a,b=chromatogram(reference(),norm);ax.plot(xr,a+b,color='#8b8581',ls='--',label='Saved total')
        ax.set(xlabel='Normalized time t/t₀' if norm else 'Time (min)',ylabel='Area-normalized signal',ylim=(0,None))
        ax.legend(fontsize=8,frameon=False);return fig

    @render.plot(alt='Resolution versus velocity with a 1.5 target and later-peak retention time on a second panel.')
    def tradeoff():
        p=settings();fig=Figure(figsize=(12,2.7),layout='constrained');axs=fig.subplots(1,2)
        u=np.linspace(.1,5,400);ms=[calculate(Conditions(**{**asdict(p),'u':v})) for v in u]
        for ax,key,label,color in zip(axs,['Rs','t2'],['Resolution Rs','Later-peak time (min)'],['#b5122a','#c57a29']):
            scale=60 if key=='t2' else 1
            ax.plot(u,[m[key]/scale for m in ms],color=color,lw=2)
            ax.scatter([p.u],[calculate(p)[key]/scale],color='#202124',zorder=5)
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
