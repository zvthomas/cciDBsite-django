# Generated by Django 3.2.15 on 2022-10-14 21:46

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interactions', '0016_pathwaycorrelations'),
    ]

    operations = [
        migrations.AddField(
            model_name='pathway',
            name='evidences',
            field=models.TextField(default=''),
            preserve_default=False,
        ),
    ]
