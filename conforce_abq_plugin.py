"""
This is an Abaqus plugin
# TODO:
"""
from __future__ import print_function

import abaqusGui as gui
import kernelAccess as ka


class MyDialog(gui.AFXDataDialog):
    def __init__(self, owner):
        gui.AFXDataDialog.__init__(
            self,
            owner,
            title="Configurational Force Computation",
            actionButtonIds=self.APPLY | self.CANCEL
        )

        self.label = gui.FXLabel(self, "")

        # odb
        gb = gui.FXGroupBox(
            self,
            text='odb:',
            opts=gui.FRAME_GROOVE
        )
        gui.FXLabel(gb, "Compute the requested field output in this odb.")
        drop_down = gui.AFXComboBox(
            gb,
            ncols=0,
            nvis=3,
            text="select odb:",
            tgt=owner.kw_odb_name
        )
        odb_names = ka.session.odbs.keys()
        if len(odb_names) == 0:
            self.label.setText("No odb found. Open an odb file. " + self.label.getText())
        else:
            owner.kw_odb_name.setValue(odb_names[0])

        for odb_name in odb_names:
            drop_down.appendItem(odb_name)

        # method
        gb = gui.FXGroupBox(
            self,
            text='computation method:',
            opts=gui.FRAME_GROOVE
        )
        gui.FXLabel(gb, "Define how CONF_STRESS and CONF_FORCE are computed.")
        drop_down = gui.AFXComboBox(
            gb,
            ncols=0,
            nvis=2,
            text="select method:",
            tgt=owner.kw_method
        )
        owner.kw_method.setValue("motion based formulation")

        drop_down.appendItem("motion based formulation")
        drop_down.appendItem("deformation based formulation")

        # energy density expression
        gb = gui.FXGroupBox(
            self,
            text='energy density:',
            opts=gui.FRAME_GROOVE
        )
        gui.FXLabel(
            gb,
            "The energy density can be any scalar field output \n"
            "that is computed at the integration points. Multiple\n"
            "field outputs and floats can be combined using: \n"
            "+, -, *, /",
            opts=gui.JUSTIFY_LEFT)
        gui.AFXTextField(
            gb,
            ncols=30,
            labelText="define energy density:",
            tgt=owner.kw_e_expression
        )
        owner.kw_e_expression.setValue("SENER+PENER")

        # request outputs
        gb = gui.FXGroupBox(
            self,
            text='field output:',
            opts=gui.FRAME_GROOVE
        )
        gui.FXLabel(gb, "Write requested field outputs into the selected odb.")
        gui.FXLabel(gb, "request field output:")
        gui.FXCheckButton(
            gb,
            text='Deformation Gradient (DEF_GRAD)',
            tgt=owner.kw_F
        )
        gui.FXCheckButton(
            gb,
            text='First Piola-Kirchhoff Stress (FIRST_PIOLA_STRESS)',
            tgt=owner.kw_P
        )
        gui.FXCheckButton(
            gb,
            text='Configurational Stress',
            tgt=owner.kw_CS
        )
        gui.AFXTextField(
            gb,
            ncols=30,
            labelText="name of configurational stress:",
            tgt=owner.kw_CS_name
        )
        owner.kw_CS_name.setValue("CONF_STRESS")
        gui.FXCheckButton(
            gb,
            text='Configurational Force',
            tgt=owner.kw_CF
        )
        owner.kw_CF.setValue(True)
        gui.AFXTextField(
            gb,
            ncols=30,
            labelText="name of configurational force:",
            tgt=owner.kw_CF_name
        )
        owner.kw_CF_name.setValue("CONF_FORCE")


class MyForm(gui.AFXForm):
    def __init__(self, owner):
        # Construct the base class.
        gui.AFXForm.__init__(self, owner)
        self.owner = owner

        self.cmd = gui.AFXGuiCommand(self, "apply", "cf_abq_main")
        self.kw_odb_name = gui.AFXStringKeyword(self.cmd, "odb_or_odb_name", True)
        self.kw_method = gui.AFXStringKeyword(self.cmd, "method", True)
        self.kw_e_expression = gui.AFXStringKeyword(self.cmd, "e_expression", True)
        self.kw_F = gui.AFXBoolKeyword(self.cmd, "request_F", True, isRequired=True)
        self.kw_P = gui.AFXBoolKeyword(self.cmd, "request_P", True, isRequired=True)
        self.kw_CS = gui.AFXBoolKeyword(self.cmd, "request_CS", True, isRequired=True)
        self.kw_CS_name = gui.AFXStringKeyword(self.cmd, "CS_name", True)
        self.kw_CF = gui.AFXBoolKeyword(self.cmd, "request_CF", True, isRequired=True)
        self.kw_CF_name = gui.AFXStringKeyword(self.cmd, "CF_name", True)

        self.dialog = None

    def getFirstDialog(self):
        self.dialog = MyDialog(self)
        return self.dialog


toolset = gui.getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerGuiMenuButton(
    object=MyForm(toolset),
    buttonText='Conf. Force',
    kernelInitString='from conforce_abq import main as cf_abq_main'
)
