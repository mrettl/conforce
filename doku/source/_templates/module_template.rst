

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}

{% block classes %}
{% if classes %}
.. rubric:: Classes

.. autosummary::
  :toctree: {{ objname }}
  :template: class2.rst
{% for item in classes %}
  {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block functions %}
{% if functions %}
.. rubric:: Functions

.. autosummary::
  :toctree: {{ objname }}
{% for item in functions %}
  {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}


{% block exceptions %}
{% if exceptions %}
.. rubric:: Exceptions

.. autosummary::
{% for item in exceptions %}
  {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
