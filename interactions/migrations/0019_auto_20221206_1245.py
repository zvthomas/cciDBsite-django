# Generated by Django 3.2.15 on 2022-12-06 20:45

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interactions', '0018_auto_20221018_0939'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pathwayandcelltype',
            name='sorr',
            field=models.CharField(choices=[('s', 'Sending'), ('r', 'Receiving')], default='s', max_length=1),
        ),
        migrations.AlterField(
            model_name='pathwaycorrelations',
            name='p1a',
            field=models.CharField(choices=[('s', 'Sending'), ('r', 'Receiving')], default='s', max_length=1),
        ),
        migrations.AlterField(
            model_name='pathwaycorrelations',
            name='p2a',
            field=models.CharField(choices=[('s', 'Sending'), ('r', 'Receiving')], default='s', max_length=1),
        ),
    ]