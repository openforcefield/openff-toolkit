{#
If parent object's name starts with a capital, include the parent
name in this name; otherwise, just this one.

This is so that methods, attributes etc get the name of their class,
but free functions do not.
#}
{%- if fullname.count(".") > 0 -%}
    {%- set namewithparent = ".".join(fullname.split(".")[-2:]) -%}
{%- endif -%}
{%- if namewithparent is defined and namewithparent[0].isupper() -%}
    {%- set title = namewithparent -%}
{%- else -%}
    {%- set title = objname -%}
{%- endif -%}
{{ ("``" ~ title ~ "``") | underline('=')}}

.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}
