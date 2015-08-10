{% macro before_run() %}
PyObject_CallMethod(_owner, "initialise_queue", NULL);
{% endmacro %}

{% macro main() %}
{% endmacro %}

{% macro support_code() %}
{% endmacro %}

{% macro python_before_main() %}
_owner.push_spikes()
{% endmacro %}
