{{ objname | escape | underline }}

.. currentmodule:: {{ module }}


.. autoclass:: {{ objname }}
	:no-members:
	
{% block methods %}

{% if methods %}

Methods
-------

.. autosummary::
	:toctree: {{ objname }}

{% for item in methods %}
	~{{ name }}.{{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block attributes %}
{% if attributes %}

Attributes
----------

.. autosummary::
	:toctree: {{ objname }}
{% for item in attributes %}
	~{{ name }}.{{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
