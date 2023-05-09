import abaqusGui as gui

toolset = gui.getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerKernelMenuButton(
    buttonText='Confor',
    moduleName='cf_abq.main', functionName='main()'
)
