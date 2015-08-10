{% macro before_run() %}
def main(ns):
    ns['_owner'].initialise_queue()
{% endmacro %}

{% macro main() %}
def main(ns):
    ns['_owner'].push_spikes()
{% endmacro %}