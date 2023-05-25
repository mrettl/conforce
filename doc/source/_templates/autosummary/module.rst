{{ objname | escape | underline }}

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}

{% block functions %}
{% if functions %}

Functions
---------

.. autosummary::
	:toctree: {{ objname }}
{% for item in functions %}
	{{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block classes %}
{% if classes %}

Classes
-------

.. autosummary::
	:toctree: {{ objname }}
{% for item in classes %}
	{{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block exceptions %}
{% if exceptions %}

Exceptions
----------

.. autosummary::
{% for item in exceptions %}
	{{ item }}
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
	{{ item }}
{%- endfor %}
{% endif %}
{% endblock %}