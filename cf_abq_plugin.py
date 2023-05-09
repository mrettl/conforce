import abaqusGui as gui

import cf_shared


toolset = gui.getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerKernelMenuButton(
    moduleName='cf_abq.main',
    functionName='main()',
    buttonText='Confor',
    version=cf_shared.version,
    author=cf_shared.author,
    description=cf_shared.description,
    helpUrl=cf_shared.helpUrl
)
