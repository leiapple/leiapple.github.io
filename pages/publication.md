---
layout: page
title: Publication
description: A publication list with full access (for persnal useage)
keywords: Pulications
comments: true
menu: 
permalink: /publication/
---

> Publications related to machine learning interactomic potential.

<ul>
{% for pubs in site.data.pubs %}
  {% if pubs.src == 'science' %}
  <li><a href="{{ link.url }}" target="_blank">{{ link.name}}</a></li>
  {% endif %}
{% endfor %}
</ul>