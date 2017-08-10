from cycler import cycler

# see http://matplotlib.org/users/customizing.html for all options

style1 = {
    # Line styles
    'lines.linewidth': 2,
    'lines.antialiased': True,

    # Font
    'font.size': 16.0,
    'font.family': 'sans-serif',

    # Axes
    'axes.linewidth': 1.5,
    'axes.titlesize': 'x-large',
    'axes.labelsize': 'medium',
    'axes.prop_cycle': cycler('color', [
        '#1f77b4',  # blue
        '#ff7f0e',  # orange
        '#2ca02c',  # green
        '#d62728',  # red
        '#9467bd',  # purple
        '#8c564b',  # brown
        '#e377c2',  # magenta
        '#7f7f7f',  # gray
        '#bcbd22',  # yellow
        '#17becf',  # cyan
    ]),

    # Ticks
    'xtick.major.size': 6,
    'xtick.minor.size': 3,
    'xtick.major.width': 1.5,
    'xtick.minor.width': 1.5,
    'xtick.major.pad': 6,
    'xtick.minor.pad': 3,
    'xtick.labelsize': 'medium',
    'xtick.direction': 'in',
    'xtick.top': True,
    'xtick.bottom': True,
    'xtick.minor.visible': True,

    'ytick.major.size': 6,
    'ytick.minor.size': 3,
    'ytick.major.width': 1.5,
    'ytick.minor.width': 1.5,
    'ytick.major.pad': 6,
    'ytick.minor.pad': 3,
    'ytick.labelsize': 'medium',
    'ytick.direction': 'in',
    'ytick.right': True,
    'ytick.left': True,
    'ytick.minor.visible': True,

    # Legend
    'legend.fancybox': False,
    'legend.fontsize': 'medium',
    'legend.scatterpoints': 1,
    'legend.numpoints': 1,
    'legend.loc': 'best',

    # Figure
    'figure.figsize': [8.5, 9],
    'figure.titlesize': 'large',

    # Images
    'image.cmap': 'magma',
    'image.origin': 'lower',

    # Saving
    'savefig.bbox': 'tight',
    'savefig.format': 'pdf',
}

# https://matplotlib.org/users/customizing.html
#>>> plt.rcParams.keys()
#[u'_internal.classic_mode', u'agg.path.chunksize', u'animation.avconv_args', u'animation.avconv_path', u'animation.bitrate', u'animation.codec', u'animation.convert_args', u'animation.convert_path', u'animation.ffmpeg_args', u'animation.ffmpeg_path', u'animation.frame_format', u'animation.html', u'animation.mencoder_args', u'animation.mencoder_path', u'animation.writer', u'axes.autolimit_mode', u'axes.axisbelow', u'axes.edgecolor', u'axes.facecolor', u'axes.formatter.limits', u'axes.formatter.offset_threshold', u'axes.formatter.use_locale', u'axes.formatter.use_mathtext', u'axes.formatter.useoffset', u'axes.grid', u'axes.grid.axis', u'axes.grid.which', u'axes.hold', u'axes.labelcolor', u'axes.labelpad', u'axes.labelsize', u'axes.labelweight', u'axes.linewidth', u'axes.prop_cycle', u'axes.spines.bottom', u'axes.spines.left', u'axes.spines.right', u'axes.spines.top', u'axes.titlepad', u'axes.titlesize', u'axes.titleweight', u'axes.unicode_minus', u'axes.xmargin', u'axes.ymargin', u'axes3d.grid', u'backend', u'backend.qt4', u'backend.qt5', u'backend_fallback', u'boxplot.bootstrap', u'boxplot.boxprops.color', u'boxplot.boxprops.linestyle', u'boxplot.boxprops.linewidth', u'boxplot.capprops.color', u'boxplot.capprops.linestyle', u'boxplot.capprops.linewidth', u'boxplot.flierprops.color', u'boxplot.flierprops.linestyle', u'boxplot.flierprops.linewidth', u'boxplot.flierprops.marker', u'boxplot.flierprops.markeredgecolor', u'boxplot.flierprops.markerfacecolor', u'boxplot.flierprops.markersize', u'boxplot.meanline', u'boxplot.meanprops.color', u'boxplot.meanprops.linestyle', u'boxplot.meanprops.linewidth', u'boxplot.meanprops.marker', u'boxplot.meanprops.markeredgecolor', u'boxplot.meanprops.markerfacecolor', u'boxplot.meanprops.markersize', u'boxplot.medianprops.color', u'boxplot.medianprops.linestyle', u'boxplot.medianprops.linewidth', u'boxplot.notch', u'boxplot.patchartist', u'boxplot.showbox', u'boxplot.showcaps', u'boxplot.showfliers', u'boxplot.showmeans', u'boxplot.vertical', u'boxplot.whiskerprops.color', u'boxplot.whiskerprops.linestyle', u'boxplot.whiskerprops.linewidth', u'boxplot.whiskers', u'contour.corner_mask', u'contour.negative_linestyle', u'datapath', u'date.autoformatter.day', u'date.autoformatter.hour', u'date.autoformatter.microsecond', u'date.autoformatter.minute', u'date.autoformatter.month', u'date.autoformatter.second', u'date.autoformatter.year', u'docstring.hardcopy', u'errorbar.capsize', u'examples.directory', u'figure.autolayout', u'figure.dpi', u'figure.edgecolor', u'figure.facecolor', u'figure.figsize', u'figure.frameon', u'figure.max_open_warning', u'figure.subplot.bottom', u'figure.subplot.hspace', u'figure.subplot.left', u'figure.subplot.right', u'figure.subplot.top', u'figure.subplot.wspace', u'figure.titlesize', u'figure.titleweight', u'font.cursive', u'font.family', u'font.fantasy', u'font.monospace', u'font.sans-serif', u'font.serif', u'font.size', u'font.stretch', u'font.style', u'font.variant', u'font.weight', u'grid.alpha', u'grid.color', u'grid.linestyle', u'grid.linewidth', u'hatch.color', u'hatch.linewidth', u'hist.bins', u'image.aspect', u'image.cmap', u'image.composite_image', u'image.interpolation', u'image.lut', u'image.origin', u'image.resample', u'interactive', u'keymap.all_axes', u'keymap.back', u'keymap.forward', u'keymap.fullscreen', u'keymap.grid', u'keymap.home', u'keymap.pan', u'keymap.quit', u'keymap.save', u'keymap.xscale', u'keymap.yscale', u'keymap.zoom', u'legend.borderaxespad', u'legend.borderpad', u'legend.columnspacing', u'legend.edgecolor', u'legend.facecolor', u'legend.fancybox', u'legend.fontsize', u'legend.framealpha', u'legend.frameon', u'legend.handleheight', u'legend.handlelength', u'legend.handletextpad', u'legend.labelspacing', u'legend.loc', u'legend.markerscale', u'legend.numpoints', u'legend.scatterpoints', u'legend.shadow', u'lines.antialiased', u'lines.color', u'lines.dash_capstyle', u'lines.dash_joinstyle', u'lines.dashdot_pattern', u'lines.dashed_pattern', u'lines.dotted_pattern', u'lines.linestyle', u'lines.linewidth', u'lines.marker', u'lines.markeredgewidth', u'lines.markersize', u'lines.scale_dashes', u'lines.solid_capstyle', u'lines.solid_joinstyle', u'markers.fillstyle', u'mathtext.bf', u'mathtext.cal', u'mathtext.default', u'mathtext.fallback_to_cm', u'mathtext.fontset', u'mathtext.it', u'mathtext.rm', u'mathtext.sf', u'mathtext.tt', u'nbagg.transparent', u'patch.antialiased', u'patch.edgecolor', u'patch.facecolor', u'patch.force_edgecolor', u'patch.linewidth', u'path.effects', u'path.simplify', u'path.simplify_threshold', u'path.sketch', u'path.snap', u'pdf.compression', u'pdf.fonttype', u'pdf.inheritcolor', u'pdf.use14corefonts', u'pgf.debug', u'pgf.preamble', u'pgf.rcfonts', u'pgf.texsystem', u'plugins.directory', u'polaraxes.grid', u'ps.distiller.res', u'ps.fonttype', u'ps.papersize', u'ps.useafm', u'ps.usedistiller', u'savefig.bbox', u'savefig.directory', u'savefig.dpi', u'savefig.edgecolor', u'savefig.facecolor', u'savefig.format', u'savefig.frameon', u'savefig.jpeg_quality', u'savefig.orientation', u'savefig.pad_inches', u'savefig.transparent', u'scatter.marker', u'svg.fonttype', u'svg.hashsalt', u'svg.image_inline', u'text.antialiased', u'text.color', u'text.dvipnghack', u'text.hinting', u'text.hinting_factor', u'text.latex.preamble', u'text.latex.preview', u'text.latex.unicode', u'text.usetex', u'timezone', u'tk.window_focus', u'toolbar', u'verbose.fileo', u'verbose.level', u'webagg.open_in_browser', u'webagg.port', u'webagg.port_retries', u'xtick.bottom', u'xtick.color', u'xtick.direction', u'xtick.labelsize', u'xtick.major.bottom', u'xtick.major.pad', u'xtick.major.size', u'xtick.major.top', u'xtick.major.width', u'xtick.minor.bottom', u'xtick.minor.pad', u'xtick.minor.size', u'xtick.minor.top', u'xtick.minor.visible', u'xtick.minor.width', u'xtick.top', u'ytick.color', u'ytick.direction', u'ytick.labelsize', u'ytick.left', u'ytick.major.left', u'ytick.major.pad', u'ytick.major.right', u'ytick.major.size', u'ytick.major.width', u'ytick.minor.left', u'ytick.minor.pad', u'ytick.minor.right', u'ytick.minor.size', u'ytick.minor.visible', u'ytick.minor.width', u'ytick.right']
