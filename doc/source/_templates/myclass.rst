{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :toctree: .

   {% for item in methods %}
      {% if not item.startswith('__') %}
      {% if item not in inherited_members %}
      ~{{ name }}.{{ item }}
      {% endif %}
      {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}
