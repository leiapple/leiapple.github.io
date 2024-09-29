---
layout: page
title: About
description: Science, Poem, and the Sun.
keywords: Lei Zhang, 张磊
comments: true
menu: About(关于)
permalink: /about/
---

I am Lei Zhang, a Ph.D. in the field of computational material science. 


## Contact

<ul>
{% for website in site.data.social %}
<li>{{website.sitename }}：<a href="{{ website.url }}" target="_blank">@{{ website.name }}</a></li>
{% endfor %}
{% if site.url contains 'mazhuang.org' %}
<li>
微信公众号：<br />
<img style="height:238px;width:210px;border:1px solid lightgrey;" src="{{ site.url }}/assets/images/qrcode.jpg" />
</li>
{% endif %}
</ul>


## Skill Keywords

{% for skill in site.data.skills %}
### {{ skill.name }}
<div class="btn-inline">
{% for keyword in skill.keywords %}
<button class="btn btn-outline" type="button">{{ keyword }}</button>
{% endfor %}
</div>
{% endfor %}
